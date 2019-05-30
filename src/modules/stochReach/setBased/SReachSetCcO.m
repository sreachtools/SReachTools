function varargout = SReachSetCcO(method_str, sys, prob_thresh, safety_tube, ...
    options)
% Compute the stochastic reach set corresponding to the stochastic reachability 
% problem of a target tube using convex chance-constraint optimization
% =============================================================================
%
% SReachSetCcO computes the open-loop controller-based underapproximative
% stochastic reach set to the problem of stochastic reachability of a target
% tube as discussed in
%
% A. Vinod and M. Oishi, "Scalable underapproximative verification of stochastic
% LTI systems using convexity and compactness," In Proc. Hybrid Syst.: Comput. &
% Ctrl., pages 1--10, 2018. HSCC 2018
%
% A. Vinod and M. Oishi, "Stochastic reachability of a target tube: Theory and
% computation," IEEE Transactions in Automatic Control, 2018 (submitted)
% https://arxiv.org/pdf/1810.05217.pdf.
%
% =============================================================================
%
% [polytope, extra_info] = SReachSetCcO(method_str, sys, prob_thresh,...
%   safety_tube, options)
% 
% Inputs:
% -------
%   method_str  - Solution technique to be used. Must be 'chance-open'
%   sys         - System description (LtvSystem/LtiSystem object)
%   prob_thresh - Probability threshold at which the set is to be constructed
%   safety_tube - Collection of (potentially time-varying) safe sets that
%                 define the safe states (Tube object)
%   options     - Collection of user-specified options for 'chance-open'. It is
%                 a struct created using SReachSetOptions. See SReachSetOptions
%                 for details on the available options.
%
% Outputs:
% --------
%   polytope   - Underapproximative polytope of dimension sys.state_dim which
%                underapproximates the stochastic reach set
%   extra_info - A list of Matlab structs that comprises of auxillary
%                information from the set computation. It has two members
%                extra_info_wmax and extra_info_cheby.  The individual structs
%                contain the following information:
%                   1. xmax - Initial state that has the maximum reach
%                             probability to stay with the safety tube using an
%                             open-loop controller (via the method in use)
%                   2. Umax - Optimal open-loop policy ((sys.input_dim) *
%                             time_horizon)-dimensional vector 
%                             U = [u_0; u_1;...; u_N] (column vector) for xmax
%                             (via the method in use)
%                   3. xmax_reach_prob 
%                           - Maximum attainable reach probability to
%                             stay with the safety tube using an open-loop
%                             controller
%                   4. opt_theta_i 
%                           - Vector comprising of scaling factors along each
%                             user-specified direction of interest
%                   5. opt_input_vec_at_vertices 
%                           - Optimal open-loop policy ((sys.input_dim) *
%                             time_horizon)-dim.  vector U = [u_0; u_1; ...;
%                             u_N] (column vector) for each vertex of the
%                             polytope
%                   6. opt_reach_prob_i
%                           - Maximum attainable reach probability to stay with
%                             the safety tube at the vertices of the polytope
%                   7. vertices_underapprox_polytope
%                           - Vertices of the polytope
%                               xmax + opt_theta_i * options.set_of_dir_vecs
%
% Notes:
% ------
% * extra_info(1).xmax_reach_prob is the highest prob_thresh that may be given
%   while obtaining a non-trivial underapproximation
% * We compute the set by ray-shooting algorithm that guarantees an
%   underapproximation due to the compactness and convexity of the stochastic
%   reach set. 
% * SReachSetCcO has the following approaches (specified via
%   options.compute_style) for computing the origin of the rays (referred to as
%   anchor):
%   'max_safe_init' --- Choose the anchor within safe set at t=0 such that an
%                       admissible open-loop controller exists which provides
%                       maximum safety
%   'cheby'         --- Choose the anchor which is the Chebyshev center of the
%                       safe set at t=0, that also admits an open-loop 
%                       stochastic reach probability, above the prescribed
%                       probability threshold.
% * SReachSetCcO also admits a options.compute_style = 'all' to perform the
%   polytope computation from all of the above methods, and compute the convex
%   hull of the union. This approach is also guaranteed to be an
%   underapproximation, due to the convexity of the open-loop stochastic reach
%   set.  However, this approach can result in twice the number of direction
%   vectors or vertices in the underapproximative polytope.
% * See @LtiSystem/getConcatMats for more information about the
%     notation used.
% 
% =============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://sreachtools.github.io/license/
%
%

    myeps = 1e-10; % Proxy for 0. Ideally, must be eps but MPT needs slack
    validatestring(method_str,{'chance-open'});
    
    inpar = inputParser();
    inpar.addRequired('sys', @(x) validateattributes(x, ...
        {'LtiSystem','LtvSystem'}, {'nonempty'}));
    inpar.addRequired('prob_thresh', @(x) validateattributes(x, {'numeric'}, ...
        {'scalar','>=',0.5,'<=',1}));
    inpar.addRequired('safety_tube',@(x) validateattributes(x,{'Tube'}, ...
        {'nonempty'}));

    try
        inpar.parse(sys, prob_thresh, safety_tube);
    catch err
        exc = SrtInvalidArgsError.withFunctionName();
        exc = exc.addCause(err);
        throwAsCaller(exc);
    end

    % Ensure that options are provided are appropriate
    otherInputHandling(method_str, sys, options);


    % Construct the constrained initial safe set
    init_safe_set = Polyhedron('H', safety_tube(1).H, ...
                               'He',[safety_tube(1).He;
                                     options.init_safe_set_affine.He]);
    
    % Compute mean_X_zizs, cov_X_sans_input
    time_horizon = length(safety_tube)-1;
    [~, ~, G] = sys.getConcatMats(time_horizon);
    GW = G * sys.dist.concat(time_horizon);
    mean_X_zizs = GW.mean();
    cov_X_sans_input = GW.cov();
    
    % Compute \sqrt{h_i^\top * \Sigma_X_no_input * h_i}
    [concat_safety_tube_A, ~] = safety_tube.concat([1 time_horizon]+1);

    % cholesky > cov_X_sans_input = sqrt_cov_X_sans_input'*sqrt_cov_X_sans_input
    [sqrt_cov_X_sans_input, p] = chol(cov_X_sans_input);
    if p > 0
        % Non-positive definite matrix can not use Cholesky's decomposition
        % Use sqrt to obtain a symmetric non-sparse square-root matrix
        sqrt_cov_X_sans_input = sqrt(cov_X_sans_input);
    end
    % Hence, the transpose is needed
    scaled_sigma_vec = norms(concat_safety_tube_A*sqrt_cov_X_sans_input', 2,2);

    %% Step 1: Find initial state (xmax) w/ max open-loop stochastic reach prob
    if options.verbose >= 1
        disp('Computing an initial state with maximum safety probability');
    end
    xmax_soln = computeWmax(sys, options, init_safe_set, prob_thresh, ...
                    safety_tube, mean_X_zizs, scaled_sigma_vec);
    % Default values for extra_info
    extra_info_cheby = create_dummy_extra_info_empty();
    extra_info_wmax = create_dummy_extra_info_empty();

    % Check if the maximally safe initial state is above the probability
    % threshold
    if ~strcmpi(xmax_soln.cvx_status, 'Solved') ||...
            xmax_soln.reach_prob < prob_thresh
        if options.verbose >= 1
            if strcmpi(xmax_soln.cvx_status, 'Solved')
                fprintf(['Maximum reach probability: %1.4f < %1.4f (given ', ...
                    'threshold) | No polytope exists \n'], ...
                    xmax_soln.reach_prob, prob_thresh);
            else
                fprintf(['No polytope exists and non-trivial Wmax ', ...
                    'computation was not possible | CVX status: %s\n'], ...
                    xmax_soln.reach_prob, prob_thresh, xmax_soln.cvx_status);
            end
        end            
        varargout{1} = Polyhedron.emptySet(2);
        if strcmpi(xmax_soln.cvx_status, 'Solved')
            % CVX solved the problem --- update extra_info_wmax
            extra_info_wmax = create_dummy_extra_info_wmax(xmax_soln);
        else
            % CVX could not solve the problem --- reuse dummy extra_info_wmax
        end
        varargout{2} = [extra_info_wmax, extra_info_cheby];
    else
        if options.verbose >= 1
            fprintf(['Maximum reach probability: %1.2f > %1.2f (given ', ...
                'threshold) | Polytope exists \n'], ...
                xmax_soln.reach_prob, prob_thresh);
        end
        % Step 2: Check if user-provided set_of_dir_vecs are in affine hull  
        % Step 2a: Compute xmax + directions \forall directions
        n_dir_vecs = size(options.set_of_dir_vecs,2);
        check_set_of_dir_vecs = repmat(xmax_soln.xmax,1,n_dir_vecs) + ...
            options.set_of_dir_vecs;
        % Step 2b: Does it satisfy the equality constraints of init_safe_set
        xmaxPlusdirs_within_init_safe_set_eq = abs(init_safe_set.Ae * ...
            check_set_of_dir_vecs - init_safe_set.be) < myeps;
        % Step 2c: If they don't, throw error
        if ~all(all(xmaxPlusdirs_within_init_safe_set_eq))
            throw(SrtInvalidArgsError(['set_of_dir_vecs+xmax does ', ... 
                'not lie in the affine hull intersecting safe_set.']));
        end

        %% Step 3: Compute convex hull of rel. boundary points by ray-shooting        
        if any(strcmpi(options.compute_style, {'all', 'max_safe_init'}))
            % Step 3a: Obtain the polytope from max safe init
            if options.verbose >= 1
                disp('Computing the polytope via maximally safe initial state');
            end
            if nargout > 1
                [polytope_wmax, extra_info_wmax] = computePolytopeFromAnchor(...
                    xmax_soln, sys, options, init_safe_set, prob_thresh, ...
                    safety_tube, mean_X_zizs, scaled_sigma_vec);
            else
                polytope_wmax = computePolytopeFromAnchor(xmax_soln, sys, ...
                    options, init_safe_set, prob_thresh, safety_tube, ...
                    mean_X_zizs, scaled_sigma_vec);
            end
        end
        if any(strcmpi(options.compute_style, {'all', 'cheby'}))
            % Step 3b: Shoot rays from Chebyshev center of init_safe_set with
            % W_0^\ast(x_cheby) >= prob_thresh
            % Step 3b-i: Find the Chebyshev center
            if options.verbose >= 1
                disp(['Computing the Chebyshev center of initial safe set ', ...
                    'with safety probability above the prescribed threshold']);
            end
            cheby_soln = computeChebyshev(sys, options, init_safe_set, ...
                           prob_thresh, safety_tube, mean_X_zizs, ...
                           scaled_sigma_vec);
            % Step 3b-ii: Find the polytope                           
            if options.verbose >= 1
                disp('Computing the polytope via the Chebyshev center');
            end
            if nargout > 1
                [polytope_cheby, extra_info_cheby] = ...
                    computePolytopeFromAnchor(cheby_soln, sys, options, ...
                        init_safe_set, prob_thresh, safety_tube, mean_X_zizs,...
                        scaled_sigma_vec);
            else
                polytope_cheby = computePolytopeFromAnchor(cheby_soln, sys, ....
                    options,init_safe_set, prob_thresh, safety_tube, ...
                    mean_X_zizs, scaled_sigma_vec);
            end
        end
        % Step 4: Convex hull of the vertices
        switch options.compute_style
            case 'max_safe_init'
                varargout{1} = polytope_wmax;                
            case 'cheby'
                varargout{1} = polytope_cheby;                
            case 'all'
                varargout{1} = Polyhedron('V',[polytope_cheby.V;
                                               polytope_wmax.V]);                
        end
        
        % Populate varargout if extra_info is requested
        if nargout > 1
            if strcmpi(options.compute_style, 'cheby') 
                % extra_info_wmax is still []. Updated it with dummy with
                % appropriate fields filled
                varargout{2} = [create_dummy_extra_info_wmax(xmax_soln), ...
                    extra_info_cheby];
            else
                % Cheby either is either
                %  not required (compute_style = 'max_safe_init'; reuse dummmy) 
                % OR
                %  computed (compute_style = 'all'; extra_info_cheby populated)
                varargout{2} = [extra_info_wmax, extra_info_cheby];
            end
        end
    end
end

function extra_info = create_dummy_extra_info_wmax(xmax_soln)
    % Create a struct for extra_info_wmax since not computed with dummy
    % values for the opt_XXXXX and vertices_underapprox_polytope
    extra_info.xmax = xmax_soln.xmax;
    extra_info.Umax = xmax_soln.Umax;
    extra_info.xmax_reach_prob = xmax_soln.reach_prob; 
    extra_info.opt_theta_i = [];
    extra_info.opt_input_vec_at_vertices = [];
    extra_info.opt_reach_prob_i = [];
    extra_info.vertices_underapprox_polytope = [];
end

function extra_info = create_dummy_extra_info_empty()
    % Create a struct for extra_info when not computed
    extra_info.xmax = [];
    extra_info.Umax = [];
    extra_info.xmax_reach_prob = []; 
    extra_info.opt_theta_i = [];
    extra_info.opt_input_vec_at_vertices = [];
    extra_info.opt_reach_prob_i = [];
    extra_info.vertices_underapprox_polytope = [];
end

function otherInputHandling(method_str, sys, options)
    % Input handling for SReachSetCcO
    
    % Consider updating SReachSetGpO.m if any changes are made here
    
    % Ensure Gaussian-perturbed system
    validateattributes(sys.dist, {'RandomVector'}, {'nonempty'},...
        'SReachSetCcO/otherInputHandling', 'sys.dist');
    validatestring(sys.dist.type, {'Gaussian'}, {'nonempty'},...;
        'SReachSetCcO/otherInputHandling', 'sys.dist.type');
    
    % Check if prob_str and method_str are consistent        
    if ~strcmpi(options.prob_str,'term')
        throwAsCaller(SrtInvalidArgsError(['Mismatch in prob_str in the ', ...
            'options']));
    end
    if ~strcmpi(options.method_str,method_str)
        throwAsCaller(SrtInvalidArgsError(['Mismatch in method_str in the ', ...
            'options']));
    end        
    
    % Make sure the user specified set_of_dir_vecs is valid
    if size(options.set_of_dir_vecs, 1) ~= sys.state_dim
        throwAsCaller(SrtInvalidArgsError(sprintf(['set_of_dir_vecs should', ...
            ' be a collection of %d-dimensional column vectors.'], ...
            sys.state_dim)));
    end
    if size(options.set_of_dir_vecs, 2) < 2, ...
        throwAsCaller(SrtInvalidArgsError(['set_of_dir_vecs should at ', ...
            'least have two directions.']));
    end
    % Make sure the user specified init_safe_set_affine is of the correct
    % dimension
    if options.init_safe_set_affine.Dim ~= sys.state_dim &&...
            ~isempty(options.init_safe_set_affine.H)
        throwAsCaller(SrtInvalidArgsError(['init_safe_set_affine must be ', ...
            'an sys.state_dim-dimensional affine set']));
    end    
end


function [polytope, extra_info] = computePolytopeFromAnchor(x_anchor, sys, ...
    options, init_safe_set, prob_thresh, safety_tube, mean_X_zizs, ...
    scaled_sigma_vec)
    % Solve the following optimization problem
    %
    % maximize theta
    % subject to
    %       W_0(boundary_point, U) >= prob_thresh
    %       boundary_point = x_anchor.xmax + theta * dir
    %       boundary_point \in TargetTube(0)\cap Init_Safe_Set_Affine
    %       U \in U^N
    %
    % The implementation of the constraint W_0(boundary_point, U) >= prob_thresh
    % is enforced by Boole's inequality. Let mu_X(U) be the mean trajectory
    % of the concatenated state and sigma_X be its covariance. Let the target 
    % tube be the polytope {a_i^T * X <= b_i} for a collection of hyperplanes.
    %
    % Then, for scaled_sigma_vec = sqrt(a_i^T * sigma_X * a_i), we have
    % W_0(boundary_point, U)>= 1 - sum(delta_i) provided
    %
    %      norminv >= c_j * delta_i + d_j                             (A)
    %      mu_X(U) = Z boundary_point + H U + G mu_W
    %      a_i^T * mu_X(U) + scaled_sigma_vec * norminv <= b_i
    %      lb_delta <= delta_i 
    %      delta_i <= prob_thresh
    %      1 - sum(delta_i) >= prob_thresh
    % 
    % Here, 
    %
    % * (A) is a convex implementation of norminv = max_j (c_j * delta_i + d_j)
    % * max_j (c_j * y + d_j) is the piecewise affine approximation of
    %   \Phi^{-1}(1-y) for y in [lb_delta,0.5] where \Phi^{-1}(.) is the 
    %   standard normal variable's inverse cumulative density function function. 
    %   This function is known to be convex, and computeNormCdfInvOverApprox 
    %   provides c_j,d_j.
    % * The definition of x_anchor is given by computeWmax or computeChebyshev.
    
    %% Repetition
    % Get time_horizon, concat_XXX_A, concat_XXX_b, n_lin_const, Z, H, 
    % invcdf_approx_{m,c}, lb_deltai
    % This entire section evaluation took only 0.05 seconds on a 2.10 GHz
    % i7 laptop
    time_horizon = length(safety_tube)-1;
    [concat_safety_tube_A, concat_safety_tube_b] =...
        safety_tube.concat([1 time_horizon]+1);
    n_lin_const = size(concat_safety_tube_A,1);
    [concat_input_space_A, concat_input_space_b] = getConcatInputSpace(sys, ...
        time_horizon);
    [Z, H, ~] = getConcatMats(sys, time_horizon);    
    [invcdf_approx_m, invcdf_approx_c, lb_deltai] =...
        computeNormCdfInvOverApprox(1 - prob_thresh, options.pwa_accuracy, ...
            n_lin_const);
    
    % Pre-allocation of relevant outputs
    n_dir_vecs = size(options.set_of_dir_vecs,2);
    opt_input_vec_at_vertices = zeros(sys.input_dim*time_horizon,n_dir_vecs);
    opt_theta_i = zeros(1, n_dir_vecs);
    opt_reach_prob_i = zeros(1, n_dir_vecs);
    vertices_underapprox_polytope = zeros(sys.state_dim, n_dir_vecs);
    
    %% Iterate over all direction vectors + xmax
    if options.verbose >= 1
        % Vector to store the unsolved directions which is reported in the end
        unsolved_directions = [];
        fprintf('Analyzing direction :%4d/%4d', 0, n_dir_vecs);
    end
    for direction_index = 1: n_dir_vecs
        % Get direction_index-th direction in the hyperplane
        direction = options.set_of_dir_vecs(:,direction_index);
        
        if options.verbose >= 1
            fprintf('\b\b\b\b\b\b\b\b\b%4d/%4d', direction_index, n_dir_vecs);
        end
    
        %% Compute the boundary point
        cvx_begin quiet
            variable mean_X(sys.state_dim * time_horizon, 1);
            variable deltai(n_lin_const, 1);
            variable norminvover(n_lin_const, 1);
            variable theta nonnegative;
            variable boundary_point(sys.state_dim, 1);
            variable U_vector(sys.input_dim * time_horizon,1);
    
            maximize theta
            subject to
                % Chance-constraint reformulation of safety prob
                % (1) Slack variables that are bounded below by the piecewise 
                % linear approximation of \phi^{-1}(1-\delta)
                for deltai_indx = 1:n_lin_const
                    norminvover(deltai_indx) >= invcdf_approx_m .* ...
                        deltai(deltai_indx) + invcdf_approx_c; 
                end
                % (2) Trajectory
                mean_X == Z * boundary_point + H * U_vector + mean_X_zizs;
                % (3) Ono's type of reformulation of chance
                % constraints
                concat_safety_tube_A * mean_X + scaled_sigma_vec .* ...
                    norminvover <= concat_safety_tube_b;
                % (4) Lower bound on delta due to pwl's domain
                deltai >= lb_deltai;
                % (5) Upper bound on delta to make \phi^{-1}(1-\delta) concave
                %     and within the PWA approximation of CDF
                deltai <= 1 - prob_thresh;
                % (6) W_0(x,U) >= alpha
                1-sum(deltai) >= prob_thresh;
                % (7) Safe boundary point
                boundary_point == x_anchor.xmax + theta * direction;
                init_safe_set.A * boundary_point <= init_safe_set.b;
                init_safe_set.Ae * boundary_point == init_safe_set.be;
                % (8) Input constraints
                concat_input_space_A * U_vector <= concat_input_space_b;                
        cvx_end
        switch cvx_status
            case 'Solved'
                opt_theta_i(direction_index) = theta;
                opt_input_vec_at_vertices(:,direction_index) = U_vector;
                opt_reach_prob_i(direction_index) = 1 - sum(deltai);
                vertices_underapprox_polytope(:, direction_index) =...
                    boundary_point;
            otherwise
                % Compute the inequality slack
                inside_slack = max(init_safe_set.A * boundary_point -...
                    init_safe_set.b);
                % Compute the equality slack
                inplane_slack = max(init_safe_set.Ae * boundary_point -...
                    init_safe_set.be);
                if inside_slack <= eps && abs(inplane_slack) < eps
                    % All constraints satisfied
                    opt_theta_i(direction_index) = theta;
                    opt_input_vec_at_vertices(:,direction_index) = U_vector;
                    opt_reach_prob_i(direction_index) = 1 - sum(deltai);
                    vertices_underapprox_polytope(:, direction_index) =...
                        boundary_point;
                else
                    warning('SReachTools:runtime', sprintf(['CVX could not ',...
                        'solve the line search problem %d/%d. CVX status',...
                        ': %s'], direction_index, n_dir_vecs, cvx_status));
                    if options.verbose == 1
                        fprintf(['CVX violated constraint requirements: ', ...
                            'Inequality: %1.2e | Equality: %1.2e\n'], ...
                            inside_slack, inplane_slack);
                        % Storing the unsolved direction for reporting
                        unsolved_directions(end+1) = direction_index;                     
                    end
                    % Default values --- theta is zero, vertex is xmax, and
                    % inputs and reach probability as NaN
                    opt_theta_i(direction_index) = 0;
                    vertices_underapprox_polytope(:, direction_index) =...
                        x_anchor.xmax;
                    opt_input_vec_at_vertices(:,direction_index) = NaN;
                    opt_reach_prob_i(direction_index) = NaN;            
                end
        end        
    end
    if options.verbose >= 1 
        fprintf('\n');
        if ~isempty(unsolved_directions)
            fprintf('Errored in %d direction vectors: ', ...
                length(unsolved_directions));
            disp(unsolved_directions);
        end
    end
    polytope = Polyhedron('V', vertices_underapprox_polytope');
    if nargout > 1
        extra_info.xmax = x_anchor.xmax;
        extra_info.Umax = x_anchor.Umax;
        extra_info.xmax_reach_prob = x_anchor.reach_prob;
        extra_info.opt_theta_i = opt_theta_i;
        extra_info.opt_input_vec_at_vertices = opt_input_vec_at_vertices;
        extra_info.opt_reach_prob_i = opt_reach_prob_i;
        extra_info.vertices_underapprox_polytope =vertices_underapprox_polytope;
    end
end


function xmax_soln = computeWmax(sys, options, init_safe_set, prob_thresh, ...
    safety_tube, mean_X_zizs, scaled_sigma_vec)
    % Solve the following optimization problem
    %
    % maximize W_0(xmax, Umax)
    % subject to
    %       xmax \in TargetTube(0) \cap Init_Safe_Set_Affine
    %       Umax \in U^N
    %
    % Using the same notation as above, we maximize (1-sum(delta_i)) with 
    % W_0(boundary_point, U)>= (1-sum(delta_i)), provided
    %
    %      norminv >= c_j * delta_i + d_j                            
    %      mu_X(U) = Z xmax + H Umax + G mu_W
    %      a_i^T * mu_X(U) + scaled_sigma_vec * norminv <= b_i
    %      lb_delta <= delta_i 
    %      delta_i <= 0.5
    %      1 - sum(delta_i) >= 0
    % 
    
    %% Repeat of above lines --- get time_horizon, concat_XXX_A, concat_XXX_b, 
    %% n_lin_const, Z, H, invcdf_approx_{m,c}, lb_deltai
    % This entire section evaluation took only 0.05 seconds on a 2.10 GHz
    % i7 laptop
    time_horizon = length(safety_tube)-1;
    [concat_safety_tube_A, concat_safety_tube_b] = ...
        safety_tube.concat([1 time_horizon]+1);
    n_lin_const = size(concat_safety_tube_A,1);
    [concat_input_space_A, concat_input_space_b] = getConcatInputSpace(sys, ...
        time_horizon);
    [Z, H, ~] = getConcatMats(sys, time_horizon);    
    [invcdf_approx_m, invcdf_approx_c, lb_deltai] =...
        computeNormCdfInvOverApprox(0.5, options.pwa_accuracy, ...
            n_lin_const);

    cvx_begin quiet
        variable mean_X(sys.state_dim * time_horizon, 1);
        variable deltai(n_lin_const, 1);
        variable norminvover(n_lin_const, 1);
        variable xmax(sys.state_dim, 1);
        variable Umax(sys.input_dim * time_horizon,1);

        maximize (1-sum(deltai))
        subject to
            % Chance-constraint reformulation of the safety prob
            % (1) Slack variables that are bounded below by the piecewise linear
            % approximation of \phi^{-1}(1-\delta)
            for deltai_indx = 1:n_lin_const
                norminvover(deltai_indx) >= invcdf_approx_m.* ...
                    deltai(deltai_indx) + invcdf_approx_c; 
            end
            % (2) Trajectory
            mean_X == Z * xmax + H * Umax + mean_X_zizs;
            % (3) Ono's type of reformulation of chance constraints
            concat_safety_tube_A* mean_X + scaled_sigma_vec .* norminvover...  
                    <= concat_safety_tube_b;
            % (4) Lower bound on delta due to pwl's domain
            deltai >= lb_deltai;
            % (5) Upper bound on delta to make \phi^{-1}(1-\delta) concave
            %     and within the PWA approximation of CDF
            deltai <= 0.5;
            % (6) W_0(x,U) >= 0
            1-sum(deltai) >= 0;
            % (7) Safety constraints on the initial state
            init_safe_set.A * xmax  <= init_safe_set.b
            init_safe_set.Ae * xmax == init_safe_set.be
            % (8) Input constraints
            concat_input_space_A * Umax <= concat_input_space_b;                        
    cvx_end
    if ~strcmpi(cvx_status, 'Solved') && ~strcmpi(cvx_status, 'Infeasible') 
        warning('SReachTools:runtime', ['CVX failed to obtain the initial',...
            ' state with maximum reach probability, potentially due to ',...
            'numerical issues.']); 
    end
    xmax_soln.reach_prob = 1-sum(deltai);
    xmax_soln.xmax = xmax;
    xmax_soln.Umax = Umax; 
    xmax_soln.cvx_status = cvx_status;
end

function cheby_soln = computeChebyshev(sys, options, init_safe_set, ...
    prob_thresh, safety_tube, mean_X_zizs, scaled_sigma_vec)
    % Solve the following optimization problem
    %
    % maximize R
    % subject to
    %       W_0(xmax, Umax) >= alpha
    %       (xmax + d)\in TargetTube(0) for every n-dim. vector d, ||d|| <= R
    %       xmax\in Init_Safe_Set_Affine 
    %       Umax \in U^N
    %
    % Note that second and third constraints together imply 
    %       xmax \in TargetTube(0) \cap Init_Safe_Set_Affine
    %
    % Using the same notation as above, we maximize (1-sum(delta_i)) with 
    % W_0(boundary_point, U)>= (1-sum(delta_i)), provided
    %
    %      norminv >= c_j * delta_i + d_j                            
    %      mu_X(U) = Z xmax + H Umax + G mu_W
    %      a_i^T * mu_X(U) + scaled_sigma_vec * norminv <= b_i
    %      lb_delta <= delta_i 
    %      delta_i <= 0.5
    %      1 - sum(delta_i) >= alpha
    % 
    
    %% Repeat of above lines --- get time_horizon, concat_XXX_A, concat_XXX_b, 
    %% n_lin_const, Z, H, invcdf_approx_{m,c}, lb_deltai
    % This entire section evaluation took only 0.05 seconds on a 2.10 GHz
    % i7 laptop
    time_horizon = length(safety_tube)-1;
    [concat_safety_tube_A, concat_safety_tube_b] =...
        safety_tube.concat([1 time_horizon]+1);
    n_lin_const = size(concat_safety_tube_A,1);
    [concat_input_space_A, concat_input_space_b] = getConcatInputSpace(sys, ...
        time_horizon);
    [Z, H, ~] = getConcatMats(sys, time_horizon);    
    [invcdf_approx_m, invcdf_approx_c, lb_deltai] =...
        computeNormCdfInvOverApprox(1-prob_thresh, options.pwa_accuracy, ...
            n_lin_const);

    % Dual norm of a_i in init_safe_set for chebyshev-centering constraint
    dual_norm_of_init_safe_set_A = norms(init_safe_set.A,2,2); 
    
    %% Compute the Chebyshev center of the W_0(x,U)>= prob_thresh
    cvx_begin quiet
        variable mean_X(sys.state_dim * time_horizon, 1);
        variable deltai(n_lin_const, 1);
        variable norminvover(n_lin_const, 1);
        variable xmax(sys.state_dim, 1);
        variable Umax(sys.input_dim * time_horizon,1);
        variable R nonnegative;
    
        maximize R
        subject to
            % Chance-constraint reformulation of safety prob
            % (1) Slack variables that are bounded below by the
            % piecewise linear approximation of \phi^{-1}(1-\delta)
            for deltai_indx = 1:n_lin_const
                norminvover(deltai_indx) >= invcdf_approx_m.* ...
                    deltai(deltai_indx) + invcdf_approx_c; 
            end
            % (2) Trajectory
            mean_X == Z * xmax + H * Umax + mean_X_zizs;
            % (3) Ono's type of reformulation of chance constraints
            concat_safety_tube_A* mean_X + scaled_sigma_vec .* norminvover...
                    <= concat_safety_tube_b;
            % (4) Lower bound on delta due to pwl's domain
            deltai >= lb_deltai;
            % (5) Upper bound on delta to make \phi^{-1}(1-\delta) concave
            %     and within the PWA approximation of CDF
            deltai <= 1 - prob_thresh;
            % (6) W_0(x,U) >= prob_thresh
            1-sum(deltai) >= prob_thresh
            % (7) Chebyshev-centering constraints-based safety
            init_safe_set.A * xmax + R * dual_norm_of_init_safe_set_A...
                        <= init_safe_set.b
            init_safe_set.Ae * xmax == init_safe_set.be
            % (8) Input constraints
            concat_input_space_A * Umax <= concat_input_space_b;
    cvx_end
    switch cvx_status
        case {'Solved','Inaccurate/Solved'}
            cheby_soln.reach_prob = 1-sum(deltai);
            cheby_soln.xmax = xmax;
            cheby_soln.Umax = Umax; 
            cheby_soln.R = R;
            cheby_soln.cvx_status = cvx_status;
        otherwise
            throw(SrtInvalidArgsError('CVX failed to obtain the ',...
                'Chebyshev center, potentially due to numerical issues.'));
    end
end
