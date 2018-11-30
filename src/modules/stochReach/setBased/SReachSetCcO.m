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
%   options     - Collection of user-specified options for 'chance-open'
%                 (Matlab struct created using SReachSetOptions)
%
% Outputs:
% --------
%   polytope   - Underapproximative polytope of dimension sys.state_dim which
%                underapproximates the stochastic reach set
%   extra_info - A list of Matlab structs that comprises of auxillary
%                information from the set computation.
%                The list has two members extra_info_wmax and extra_info_cheby.
%                Both of these structs contain the following information:
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
% * See @LtiSystem/getConcatMats for more information about the
%     notation used.
% 
% =============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
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
    mean_X_zizs = GW.parameters.mean;
    cov_X_sans_input = GW.parameters.covariance;
    
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
    sigma_vec = norms(concat_safety_tube_A*sqrt_cov_X_sans_input', 2,2);

    %% Step 1: Find initial state (xmax) w/ max open-loop stochastic reach prob
    xmax_soln = computeWmax(sys, options, init_safe_set, prob_thresh, ...
                    safety_tube, mean_X_zizs, sigma_vec);
    if ~strcmpi(xmax_soln.cvx_status, 'Solved') ||...
            xmax_soln.reach_prob < prob_thresh
        polytope = Polyhedron.emptySet(2);
        if nargout > 1
            extra_info_wmax.xmax = xmax_soln.xmax;
            extra_info_wmax.Umax = xmax_soln.Umax;
            extra_info_wmax.xmax_reach_prob = xmax_soln.reach_prob;
            extra_info_cheby = [];
        end
    else
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
        % Step 3a: Shoot rays from xmax
        if options.verbose >= 1
            disp('Computing the polytope via a maximally safe initial state');
        end
        if nargout > 1
            [polytope_wmax, extra_info_wmax] = computePolytopeFromXmax(...
                xmax_soln, sys, options, init_safe_set, prob_thresh, ...
                safety_tube, mean_X_zizs, sigma_vec);
        else
            polytope_wmax = computePolytopeFromXmax(xmax_soln, sys, options, ...
                init_safe_set, prob_thresh, safety_tube, mean_X_zizs, ...
                sigma_vec);
        end
        % Step 3b: Shoot rays from Chebyshev center of init_safe_set with
        % W_0^\ast(x_cheby) >= prob_thresh
        % Step 3b-i: Find the Chebyshev center
        cheby_soln = computeChebyshev(sys, options, init_safe_set, ...
                       prob_thresh, safety_tube, mean_X_zizs, sigma_vec);
        % Step 3b-ii: Find the polytope                           
        if options.verbose >= 1
            disp('Computing the polytope via the Chebyshev center');
        end
        if nargout > 1
            [polytope_cheby, extra_info_cheby] = computePolytopeFromXmax( ...
                cheby_soln, sys, options, init_safe_set, prob_thresh, ...
                safety_tube, mean_X_zizs, sigma_vec);
        else
            polytope_cheby = computePolytopeFromXmax(cheby_soln, sys, ....
                options,init_safe_set, prob_thresh, safety_tube, ...
                mean_X_zizs, sigma_vec);
        end
        % Step 4: Convex hull of the vertices
        polytope = Polyhedron('V',[polytope_cheby.V;
                                   polytope_wmax.V]);
    end
    varargout{1} = polytope;
    if nargout > 1
        varargout{2} = [extra_info_wmax, extra_info_cheby];
    end
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
        throwAsCaller(SrtInvalidArgsError(['set_of_dir_vecs should be a ', ...
            'collection of n-dimensional column vectors.']));
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


function [polytope, extra_info] = computePolytopeFromXmax(xmax_soln, sys, ...
    options, init_safe_set, prob_thresh, safety_tube, mean_X_zizs, ...
    sigma_vec)

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
    
    % Pre-allocation of relevant outputs
    n_dir_vecs = size(options.set_of_dir_vecs,2);
    opt_input_vec_at_vertices = zeros(sys.input_dim*time_horizon,n_dir_vecs);
    opt_theta_i = zeros(1, n_dir_vecs);
    opt_reach_prob_i = zeros(1, n_dir_vecs);
    vertices_underapprox_polytope = zeros(sys.state_dim, n_dir_vecs);
    
    %% Iterate over all direction vectors + xmax
    if options.verbose >= 1
        unsolved_directions = [];
    end
    for direction_index = 1: n_dir_vecs
        % Get direction_index-th direction in the hyperplane
        direction = options.set_of_dir_vecs(:,direction_index);
        
        if options.verbose >= 1
            fprintf('Analyzing direction :%2d/%2d\n', direction_index, ...
                n_dir_vecs);
        end
    
        %% Solve the optimization problem to compute the boundary
        %% point of the underapproximative stochastic reach set
        % maximize theta
        % subject to
        %   boundary_point = xmax + theta * direction
        %   U\in U^N
        %   boundary_point\in InitState
        %   W_0( boundary_point, U) >= prob_thresh 
        %              [As before, the computation of W_0(x,U) is
        %              done via chance-constraint reformulation and
        %              risk allocation]
        cvx_begin quiet
            variable mean_X(sys.state_dim * time_horizon, 1);
            variable deltai(n_lin_const, 1);
            variable norminvover(n_lin_const, 1);
            variable theta nonnegative;
            variable boundary_point(sys.state_dim, 1);
            variable U_vector(sys.input_dim * time_horizon,1);
    
            maximize theta
            subject to
                % Define boundary point
                boundary_point == xmax_soln.xmax + theta * direction;
                % Safe boundary point
                init_safe_set.A * boundary_point <= init_safe_set.b;
                init_safe_set.Ae * boundary_point == init_safe_set.be;
                % Input constraints
                concat_input_space_A * U_vector <= concat_input_space_b;
                % Chance-constraint reformulation of safety
                % prob(1-6) (1) Slack variables that are bounded
                % below by the piecewise linear approximation of
                % \phi^{-1}(1-\delta)
                for deltai_indx = 1:n_lin_const
                    norminvover(deltai_indx) >= invcdf_approx_m .* ...
                        deltai(deltai_indx) + invcdf_approx_c; 
                end
                % (2) Trajectory
                mean_X == Z * boundary_point + H * U_vector + mean_X_zizs;
                % (3) Ono's type of reformulation of chance
                % constraints
                concat_safety_tube_A * mean_X + sigma_vec .* ...
                    norminvover <= concat_safety_tube_b;
                % (4) Lower bound on delta due to pwl's domain
                deltai >= lb_deltai;
                % (5) Upper bound on delta to make \phi^{-1}(1-\delta) concave
                deltai <= 0.5;
                % (6) W_0(x,U) >= alpha
                1-sum(deltai) >= prob_thresh
        cvx_end
        if strcmpi(cvx_status,'Solved')
            opt_theta_i(direction_index) = theta;
            opt_input_vec_at_vertices(:,direction_index) = U_vector;
            opt_reach_prob_i(direction_index) = 1 - sum(deltai);
            vertices_underapprox_polytope(:, direction_index)= boundary_point;
        else
            warning('SReachTools:runtime', sprintf(['CVX could not solve ',...
                'the line search problem %d/%d. CVX status: %s'],...
                direction_index, n_dir_vecs, cvx_status));
            if options.verbose == 1
                unsolved_directions(end+1) = direction_index; 
                inside_slack = max(init_safe_set.A * boundary_point -...
                    init_safe_set.b);
                inplane_slack = max(init_safe_set.Ae * boundary_point -...
                    init_safe_set.be);
                if inside_slack <= eps && abs(inplane_slack) < eps
                    % All constraints satisfied
                else
                    fprintf(['CVX violated constraint requirements: Inequality:',...
                        ' %1.2e | Equality: %1.2e\n'], inside_slack, inplane_slack);
                end
            end
            % Default values --- theta is zero, vertex is xmax, and inputs and
            % reach probability as NaN
            opt_theta_i(direction_index) = 0;
            vertices_underapprox_polytope(:, direction_index)= xmax_soln.xmax;
            opt_input_vec_at_vertices(:,direction_index) = NaN;
            opt_reach_prob_i(direction_index) = NaN;            
        end        
    end
    if options.verbose >= 1 && length(unsolved_directions) > 0
        fprintf('Errored in %d direction vectors:',length(unsolved_directions));
        disp(unsolved_directions);
    end
    polytope = Polyhedron('V', vertices_underapprox_polytope');
    if nargout > 1
        extra_info.xmax = xmax_soln.xmax;
        extra_info.Umax = xmax_soln.Umax;
        extra_info.xmax_reach_prob = xmax_soln.reach_prob;
        extra_info.opt_theta_i = opt_theta_i;
        extra_info.opt_input_vec_at_vertices = opt_input_vec_at_vertices;
        extra_info.opt_reach_prob_i = opt_reach_prob_i;
        extra_info.vertices_underapprox_polytope =vertices_underapprox_polytope;
    end
end


function xmax_soln = computeWmax(sys, options, init_safe_set, prob_thresh, ...
    safety_tube, mean_X_zizs, sigma_vec)

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
        computeNormCdfInvOverApprox(1-prob_thresh, options.pwa_accuracy, ...
            n_lin_const);


    % Maximum value for the open-loop-based control
    % maximize W_0(x,U) 
    % subject to
    %   W_0(x,U) >= prob_thresh  [Computation of W_0(x,U) is done via
    %                             chance-constraint reformulation and risk
    %                             allocation]
    %   x\in InitState
    %   U\in U^N
    cvx_begin quiet
        variable mean_X(sys.state_dim * time_horizon, 1);
        variable deltai(n_lin_const, 1);
        variable norminvover(n_lin_const, 1);
        variable xmax(sys.state_dim, 1);
        variable Umax(sys.input_dim * time_horizon,1);

        maximize (1-sum(deltai))
        subject to
            % Chance-constraint reformulation of the safety prob (1-6)
            % (1) Slack variables that are bounded below by the piecewise linear
            % approximation of \phi^{-1}(1-\delta)
            for deltai_indx = 1:n_lin_const
                norminvover(deltai_indx) >= invcdf_approx_m.* ...
                    deltai(deltai_indx) + invcdf_approx_c; 
            end
            % (2) Trajectory
            mean_X == Z * xmax + H * Umax + mean_X_zizs;
            % (3) Ono's type of reformulation of chance constraints
            concat_safety_tube_A* mean_X + sigma_vec .* norminvover...  
                    <= concat_safety_tube_b;
            % (4) Lower bound on delta due to pwl's domain
            deltai >= lb_deltai;
                % (5) Upper bound on delta to make \phi^{-1}(1-\delta) concave
            deltai <= 0.5;
            % (6) Safety constraints on the initial state
            init_safe_set.A * xmax  <= init_safe_set.b
            init_safe_set.Ae * xmax == init_safe_set.be
            % (7) Input constraints
            concat_input_space_A * Umax <= concat_input_space_b;
    cvx_end
    if ~strcmpi(cvx_status, 'Solved')
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
    prob_thresh, safety_tube, mean_X_zizs, sigma_vec)
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

    %% Compute the Chebyshev center of the W_0(x,U)>= prob_thresh
    % maximize R 
    % subject to
    %   W_0(x,U) >= W_0(x_max,U_max)
    %   x\in InitState and and admits a radius R ball fit inside the
    %       inequalities
    %   U\in U^N
    % Dual norm of a_i in init_safe_set for chebyshev-centering
    % constraint
    dual_norm_of_init_safe_set_A = norms(init_safe_set.A,2,2); 
    cvx_begin quiet
        variable mean_X(sys.state_dim * time_horizon, 1);
        variable deltai(n_lin_const, 1);
        variable norminvover(n_lin_const, 1);
        variable xmax(sys.state_dim, 1);
        variable Umax(sys.input_dim * time_horizon,1);
        variable R nonnegative;
    
        maximize R
        subject to
            % Chance-constraint reformulation of safety prob (1-6)
            % (1) Slack variables that are bounded below by the
            % piecewise linear approximation of \phi^{-1}(1-\delta)
            for deltai_indx = 1:n_lin_const
                norminvover(deltai_indx) >= invcdf_approx_m.* ...
                    deltai(deltai_indx) + invcdf_approx_c; 
            end
            % (2) Trajectory
            mean_X == Z * xmax + H * Umax + mean_X_zizs;
            % (3) Ono's type of reformulation of chance constraints
            concat_safety_tube_A* mean_X + sigma_vec .* norminvover...
                    <= concat_safety_tube_b;
            % (4) Lower bound on delta due to pwl's domain
            deltai >= lb_deltai;
            % (5) Upper bound on delta to ensure \phi^{-1}(1-\delta) is
            % concave
            deltai <= 0.5;
            % (6) W_0(x,U) >= W_0(x_max,U_max)  %TODO: Or should we use Wmax?
            1-sum(deltai) >= prob_thresh
            % Chebyshev-centering constraints-based safety
            init_safe_set.A * xmax + R * dual_norm_of_init_safe_set_A...
                        <= init_safe_set.b
            init_safe_set.Ae * xmax == init_safe_set.be
            % Input constraints
            concat_input_space_A * Umax <= concat_input_space_b;
    cvx_end
    if ~strcmpi(cvx_status, 'Solved')
        warning('SReachTools:runtime', ['CVX failed to obtain the Chebyshev',...
            ' centered trajectory, potentially due to numerical issues.']); 
    end
    cheby_soln.reach_prob = 1-sum(deltai);
    cheby_soln.xmax = xmax;
    cheby_soln.Umax = Umax; 
    cheby_soln.R = R;
    cheby_soln.cvx_status = cvx_status;
end
