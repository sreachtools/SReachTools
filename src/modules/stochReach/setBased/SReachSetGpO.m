function varargout = SReachSetGpO(method_str, sys, prob_thresh, safety_tube, ...
    options)
% Compute the stochastic reach set corresponding to the stochastic reachability 
% problem of a target tube using Genz's algorithm and MATLAB's patternsearch
% =============================================================================
%
% SReachSetGpO computes the open-loop controller-based underapproximative
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
% [polytope, extra_info] = SReachSetGpO(method_str, sys, prob_thresh,...
%   safety_tube, options)
% 
% Inputs:
% -------
%   method_str  - Solution technique to be used. Must be 'genzps-open'
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
%   extra_info - A Matlab struct that comprises of auxillary
%                information from the set computation:
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
%                   8. polytope_cc_open
%                           - Polytope obtained by running SReachSet with
%                             'chance-open'
%                   9. extra_info_cco
%                           - Extra information obtained by running SReachSet 
%                             with 'chance-open'. See SReachSetCcO for more
%                             details.
% Notes:
% ------
% * extra_info.xmax_reach_prob is the highest prob_thresh that may be given
%   while obtaining a non-trivial underapproximation
% * This function calls SReachSet (chance-open) for the purposes of
%   initialization of the nonlinear solver as well as constructing bounds on the
%   line search.
% * We compute the set by ray-shooting algorithm that guarantees an
%   underapproximation due to the compactness and convexity of the stochastic
%   reach set. 
%   - Xmax computation is done with a heuristic of centering the initial
%     state with respect to the polytope obtained via SReachSet (chance-open),
%     followed by a patternsearch-based xmax computation
%   - Line search is done via bisection on theta, the scaling along the
%     direction vector of interest. 
%       - The underapproximative polytope returned by SReachSet (chance-open) is
%         utilized to compute a lower bound on theta
%       - The upper bound on the theta is obtained by the support vector of
%         the target tube at t=0.
%   - Internal computation is done using SReachPointGpO to maximize code
%     modularity and reuse.
%   - The bisection is guided by feasibility of finding an open-loop controller
%     at each value of theta that minimizes 
%
%         -log(min(Reach_prob, prob_thresh)) (1)
%
%     Reach_prob is a log-concave function over the open-loop and thus (1) is a
%     convex objective. 
%       1. Using options.thresh, an inner min operation is used that is
%          convexity-preserving. This is an attempt to ensure that the
%          Monte-Carlo simulation-driven optimization does not spend too much 
%          time looking for global optimality, when a certificate of exceeding a
%          lower bound suffices. 
%       2. After the optimization, the optimal value is reevaluated using a
%          fresh set of particles for generality.
%       3. At each value of theta, feasibility is determined if the optimal
%          value of the optimization is above prob_thresh. However, by using a
%          sufficiently low desired_accuracy in Genz's algorithm, it can be seen
%          that we will not be too off. Moreover, Hoeffding's inequality can be
%          utilized to bound the overapproximation.
%     The first two adjustments are implemented in SReachPointGpO, the
%     point-based stochastic reachability computation using genzps-open
%     method. 
% * Xmax computation is done with a heuristic of centering the initial
%   state with respect to the polytope obtained via SReachSet
%   (chance-open), followed by a patternsearch-based xmax computation
% * In 'genzps-open', desired accuracy is the farthest lower bound on the
%   confidence interval acceptable. In order to remain conservative,
%   RandomVector/getProbPolyhedron subtracts desired_accuracy from the result to
%   yield an underapproximation. For higher desired_accuracy, the result may be
%   more conservative but faster. For lower desired_accuracy, the result may
%   take more time.
% * See @LtiSystem/getConcatMats for more information about the
%     notation used.
%
% =============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://sreachtools.github.io/license/
%

    validatestring(method_str,{'genzps-open'});
    
    inpar = inputParser();
    inpar.addRequired('sys', @(x) validateattributes(x, ...
        {'LtiSystem','LtvSystem'}, {'nonempty'}));
    inpar.addRequired('prob_thresh', @(x) validateattributes(x, {'numeric'}, ...
        {'scalar','>=',0,'<=',1}));
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
    
    % Target tubes has polyhedra T_0, T_1, ..., T_{time_horizon}
    time_horizon = length(safety_tube)-1;
    
    if options.verbose >= 1 
        disp(['Compute an underapproximation of the set using chance ',...
            'constraints (SReachSetCcO)']);
    end
    cc_options = SReachSetOptions('term', 'chance-open', 'verbose', ...
        options.verbose, 'set_of_dir_vecs', options.set_of_dir_vecs, ...
        'init_safe_set_affine', options.init_safe_set_affine, ...
        'compute_style', options.compute_style_ccc);
    if nargout > 1
        % Require extra_info
        [polytope_cc_open, extra_info_cco] = SReachSet('term','chance-open', ...
            sys, prob_thresh, safety_tube, cc_options);  
    else
        % Skip extra_info
        polytope_cc_open = SReachSet('term','chance-open', sys, prob_thresh, ...
            safety_tube, cc_options);  
    end
    
    % Construct the constrained initial safe set
    init_safe_set = Polyhedron('H', safety_tube(1).H, ...
                      'He',[safety_tube(1).He;options.init_safe_set_affine.He]);
    
    if options.verbose >= 1 
        disp('Compute an initial guess for the xmax');
    end
    % Compute the chebyshev center of the initial safe set and seek the
    % optimal xmax that maximizes the open-loop safety from there
    [xmax_reach_prob, Umax, xmax] = computeXmaxViaPatternsearch( ...
        init_safe_set, extra_info_cco, sys, time_horizon, safety_tube, ... 
        options);
    
    % No. of directions
    n_dir_vecs = size(options.set_of_dir_vecs, 2);
    
    % If non-trivial underapproximative stochastic reach polytope
    if xmax_reach_prob < prob_thresh
        % Stochastic reach underapproximation is empty and no
        % admissible open-loop policy exists
        if options.verbose >= 1
            fprintf(['Polytopic underapproximation does not exist for ', ...
                 'alpha = %1.2f since W(x_max) = %1.3f.\n\n'], ...
                 prob_thresh, xmax_reach_prob);
        end

        % Assigning the outputs to trivial results
        polytope = Polyhedron();
        opt_input_vec_at_vertices =nan(sys.input_dim * time_horizon,n_dir_vecs);
        opt_theta_i = zeros(1, n_dir_vecs);
        opt_reach_prob_i = zeros(1, n_dir_vecs);
        vertices_underapprox_polytope = nan(sys.state_dim, n_dir_vecs);
    else
        if options.verbose >= 1
            % Stochastic reach underapproximation is non-trivial
            fprintf(['Polytopic underapproximation exists for alpha = ', ...
                     '%1.2f since W(x_max) = %1.3f.\n'], ...
                     prob_thresh, xmax_reach_prob);
        end

        % For storing boundary points
        opt_theta_i = zeros(1, n_dir_vecs);
        opt_reach_prob_i = zeros(1, n_dir_vecs);
        opt_input_vec_at_vertices =nan(sys.input_dim * time_horizon,n_dir_vecs);
        optionsGpO = SReachPointOptions('term', 'genzps-open', ...
            'thresh', prob_thresh, 'PSOptions', options.PSoptions, ...
            'desired_accuracy', options.desired_accuracy);
        
        %% Iterate over all direction vectors + xmax and perform a line search 
        %% via bisection
        for dir_index = 1: n_dir_vecs
            % Get direction_index-th direction in the hyperplane
            dir_vec = options.set_of_dir_vecs(:,dir_index);
            
            if options.verbose >= 1
                fprintf('\nAnalyzing direction: %2d/%2d\n',dir_index, ...
                    n_dir_vecs);
            end

            %% Bounds on theta \in [lb_theta, ub_theta] + other bisection params
            [lb_theta, ub_theta, opt_theta_so_far, opt_inputs_so_far,...
                opt_reachAvoid_prob_so_far] = computeBisectionParms(dir_vec,...
                    init_safe_set, polytope_cc_open, xmax, sys, safety_tube, ...
                    options);
            if opt_reachAvoid_prob_so_far < 0
                % Chance-constrained based initialization failed
                % Should happen only if lb_theta = 0 TODO: Check this
                lb_theta = 0;
                opt_theta_so_far = 0;
                opt_inputs_so_far = Umax;
                opt_reachAvoid_prob_so_far = xmax_reach_prob;
            end
            
            if options.verbose >= 1  && ((ub_theta - lb_theta) > ... 
                    options.tol_bisect)
                %10 chars b/n | |
                fprintf(['OptRAProb |  OptTheta  |  LB_theta  |',...
                          '  UB_theta  |  OptInp^2  | Exit reason\n']); 
            end
                
            %% Bisection-based computation of the maximum extension of
            %% the ray originating from xmax                
            % While the gap remains above tolerance_bisection
            while (ub_theta - lb_theta) > options.tol_bisect
                % Set phi as the middle point of the interval
                phi = (ub_theta+lb_theta)/2;
                % Candidate initial state provides mean_X_sans_input
                % See @LtiSystem/getConcatMats for notation
                candidate_initial_state = xmax + phi * dir_vec;
                [prob_value_at_this_step, opt_input_at_this_step] = ...
                    SReachPointGpO(sys, candidate_initial_state, safety_tube,...
                        optionsGpO);
                if prob_value_at_this_step >= prob_thresh
                    % Update values only if it is above the thresh
                    opt_reachAvoid_prob_so_far = prob_value_at_this_step;
                    opt_inputs_so_far = opt_input_at_this_step;
                    opt_theta_so_far = phi;
                    exitmessage = ' Feasible';
                else
                    % Update nothing since Genz+patternsearch failed to
                    % converge
                    exitmessage = sprintf(' Infeasible (%1.3f)', ...
                        prob_value_at_this_step);
                end
                if options.verbose >= 1
                    % Print current solution
                    fprintf(strcat('  %1.4f  |  %1.2e  |  %1.2e  |  %1.2e  ',...
                        '  |  %1.2e  |  ', exitmessage, ' \n'), ...
                        opt_reachAvoid_prob_so_far, opt_theta_so_far, ...
                        lb_theta,ub_theta,opt_inputs_so_far'*opt_inputs_so_far);
                end
                % Update bounds
                if  strcmpi(exitmessage, ' Feasible')
                    % Theta proved to be feasible
                    lb_theta = phi;
                else
                    % Shrink the upper bound
                    ub_theta = phi;
                end
            end        
            opt_theta_i(dir_index) = opt_theta_so_far;
            opt_input_vec_at_vertices(:,dir_index) = opt_inputs_so_far;
            opt_reach_prob_i(dir_index) = opt_reachAvoid_prob_so_far;
        end
        %% Construction of underapprox_stoch_reach_avoid_polytope
        % Convex hull of the boundary points 
        vertices_underapprox_polytope = xmax + ...
            opt_theta_i.*options.set_of_dir_vecs;
        polytope = Polyhedron('V', vertices_underapprox_polytope');
    end
    varargout{1} = polytope;
    if nargout > 1
        % Require extra_info
        extra_info.xmax = xmax;
        extra_info.Umax = Umax;
        extra_info.xmax_reach_prob = xmax_reach_prob;
        extra_info.opt_theta_i = opt_theta_i;
        extra_info.opt_input_vec_at_vertices = opt_input_vec_at_vertices;
        extra_info.opt_reach_prob_i = opt_reach_prob_i;
        extra_info.vertices_underapprox_polytope= vertices_underapprox_polytope;
        extra_info.polytope_cc_open = polytope_cc_open;
        extra_info.extra_info_cco = extra_info_cco;
        varargout{2} = extra_info;
    end    
end

function otherInputHandling(method_str, sys, options)
    % Input handling for SReachSetGpO
    
    % Consider updating SReachSetCcO.m if any changes are made here
    
    % Ensure Gaussian-perturbed system
    validateattributes(sys.dist, {'RandomVector'}, {'nonempty'});
    validatestring(sys.dist.type, {'Gaussian'}, {'nonempty'});
    
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
        throwAsCaller(SrtInvalidArgsError(sprintf(['init_safe_set_affine ', ...
            'must be a %d-dimensional affine set'], sys.state_dim)));
    end    
end

function [lb_theta, ub_theta, opt_theta_so_far, opt_inputs_so_far, ...
    opt_reachAvoid_prob_so_far] = computeBisectionParms(dir_vec, ...
        init_safe_set, polytope_cc_open, xmax, sys, safety_tube, options)
    %
    % Compute the bisection parameters:
    % lb_theta         - Lower bound on theta via a ray shooting on the 
    %                    chance-constrained polytope
    % ub_theta         - Upper bound on theta via a ray shooting on the safe set 
    %                    for the initial state
    % opt_theta_so_far - Same as lb_theta
    % opt_inputs_so_far, opt_reachAvoid_prob_so_far
    %                  - Obtained via chance-constrained SReachPoint at
    %                    xmax + lb_theta * dir_vec
    %
    % Notes:
    % ------
    % - The reach probability this function provides utilizes SReachPoint with
    %   chance-open for quick turnaround. For this reason, it is possible that
    %   the reach probability returned is -1 (chance-open is infeasible).
    % 


    %% Computation of the upper bound ub_theta using cvx optimization
    cvx_begin quiet
        variable theta(1) nonnegative;        
        variable safe_boundary_point(sys.state_dim, 1);
        maximize theta
        subject to
            % Define boundary point using theta
            safe_boundary_point == xmax + theta *  dir_vec;
            % Define boundary point as a convex combination of the vertices
            % defining the initial safe set (upper bound on theta)
            init_safe_set.A * safe_boundary_point <= init_safe_set.b;
    cvx_end
    switch cvx_status
        case 'Solved'
            ub_theta = theta;
        otherwise
            throw(SrtInvalidArgsError(sprintf('CVX failed (status: ',...
            '%s) to obtain the upper bound on theta for direction',...
            ' %d/%d, potentially due to numerical issues.'), ...
            cvx_status, dir_index, n_dir_vecs));
    end            
    %% Computation of the lower bound ub_theta using cvx optimization
    % Lower bound>= 0 as xmax could be a vertex. We will use
    % polytope_cc_open to further fine tune
    lb_theta = 0;
    opt_theta_so_far = 0;
    if ~polytope_cc_open.isEmptySet()
        cvx_begin quiet
            variable theta(1) nonnegative;
            variable ccc_boundary_point(sys.state_dim, 1);
            maximize theta
            subject to
                % Can't exceed the upper bound
                theta <= ub_theta; 
                % Define boundary point using theta
                ccc_boundary_point == xmax + theta * dir_vec;
                % Define boundary point as a convex combination of the vertices
                % defining polytope_cc_open (lower bound on theta)
                polytope_cc_open.A * ccc_boundary_point <= polytope_cc_open.b;
        cvx_end
        switch cvx_status
            case 'Solved'
                lb_theta = theta;
                opt_theta_so_far = lb_theta;
            otherwise
                % Both of these are set to zero by default           
                throw(SrtInvalidArgsError(sprintf('CVX failed (status: ',...
                    '%s) to obtain the lower bound on theta for direction',...
                    ' %d/%d, potentially due to numerical issues.'), ...
                    cvx_status, dir_index, n_dir_vecs));
        end
    end
    if options.verbose >= 1
        fprintf(['\b | Theta --- lower bound: %1.2f, upper bound: ',...
            '%1.2f\n'], lb_theta, ub_theta);
    end
    % Initialize the solution --- use default options
    optionsCcO = SReachPointOptions('term', 'chance-open');
    [opt_reachAvoid_prob_so_far, opt_inputs_so_far] = SReachPointCcO(...
        sys, ccc_boundary_point, safety_tube, optionsCcO);       
end

function [xmax_reach_prob, Umax, xmax] = computeXmaxViaPatternsearch(...
    init_safe_set, extra_info_ccc, sys, time_horizon, safety_tube, options)
    %
    % Compute the xmax for the patternsearch
    % 1. Compute a chebyshev center for the polytope_cc_open such that a
    %    safe mean trajectory exists
    % 2. Use the point as the initial point to compute the maximally safe
    %    initial state
    % Step 1 ensures that we are "centered".
    
    % Half-space representation for the safety tube
    [concat_safety_tube_A, concat_safety_tube_b] =...
        safety_tube.concat([1 time_horizon]+1);
    prob_poly = Polyhedron('A', concat_safety_tube_A, ...
        'b', concat_safety_tube_b);
    
    % Half-space representation for the concatenated input space
    [concat_input_space_A, concat_input_space_b] = getConcatInputSpace(sys, ...
        time_horizon);
    
    % Compute the concatenated state dynamic matrices
    [Z, H, G] = getConcatMats(sys, time_horizon); 
    
    % Compute the mean of G*W vector --- mean of concatenated state vector
    % under zero input and zero state conditions
    GW = G * sys.dist.concat(time_horizon);
    
    if options.verbose >= 1 
        fprintf(['Compute an initialization for patternsearch using ',...
            'the results from chance-open computation\n']);
    end       
    switch options.compute_style_ccc
        case 'max_safe_init'
            initial_guess_input_vector_and_xmax = [extra_info_ccc(1).Umax;
                                                   extra_info_ccc(1).xmax];
        case 'cheby'
            initial_guess_input_vector_and_xmax = [extra_info_ccc(2).Umax;
                                                   extra_info_ccc(2).xmax];
        case 'all'
            % Compute the average of the above two methods for the initial state
            % Convexity ensures that xmax guess is ok to begin with
            initial_guess_state_vector = (extra_info_ccc(1).xmax + ...
                                          extra_info_ccc(2).xmax)/2;
            initial_guess_input_vector = (extra_info_ccc(1).Umax + ...
                                          extra_info_ccc(2).Umax)/2;
            initial_guess_input_vector_and_xmax = [initial_guess_input_vector;
                                                   initial_guess_state_vector];
    end
    if isempty(initial_guess_input_vector_and_xmax)
        % In case, it is empty, set it to zero
        centered_input = Polyhedron('A', concat_input_space_A, ...
            'b', concat_input_space_b).chebyCenter();
        centered_state = init_safe_set.chebyCenter();
        initial_guess_input_vector_and_xmax = [centered_input.x;
                                               centered_state.x];
    end
    % Construct the reach cost function: 
    %   -log( max( options.desired_accuracy, ReachAvoidProbability(x_0,U)))
    % input_vector_and_xmax has first the unrolled input vector followed by
    % the initial state that is optimized | RandomVector/getProbPolyhedron
    % returns an underapproximation of the true reach probability
    negLogReachProbGivenInputVecInitState = @(input_vector_and_xmax)...
      -log( max(options.desired_accuracy, getProbPolyhedron(...
        Z * input_vector_and_xmax(sys.input_dim * time_horizon + 1: end) +...
        H * input_vector_and_xmax(1: sys.input_dim * time_horizon) + GW, ...
        prob_poly, options.desired_accuracy)));
    
    % Constraint generation --- decision variable [input_vector;xmax]
    input_vector_augmented_affine_hull_Aeq = ...
        [zeros(size(init_safe_set.Ae,1), sys.input_dim * time_horizon), ...
         init_safe_set.Ae];
    input_vector_augmented_affine_hull_beq = init_safe_set.be;
    input_vector_augmented_safe_Aineq = blkdiag(concat_input_space_A, ...
                                                init_safe_set.A);
    input_vector_augmented_safe_bineq = [concat_input_space_b;
                                         init_safe_set.b];

    % Compute xmax, the input policy, and the max reach probability
    if options.verbose >= 1 
        fprintf('Compute the xmax via patternsearch\n');
    end
    opt_input_vec_and_xmax = patternsearch(...
        negLogReachProbGivenInputVecInitState, ...
        initial_guess_input_vector_and_xmax, ...
        input_vector_augmented_safe_Aineq, ...
        input_vector_augmented_safe_bineq, ...
        input_vector_augmented_affine_hull_Aeq, ...
        input_vector_augmented_affine_hull_beq, [],[],[], options.PSoptions);
    
    %% Parse the output of patternsearch
    % Maximum attainable terminal hitting-time stochastic reach
    % probability using open-loop controller
    xmax_reach_prob = exp(-negLogReachProbGivenInputVecInitState(...
        opt_input_vec_and_xmax));
    % Optimal open_loop_control_policy
    Umax=opt_input_vec_and_xmax(1:sys.input_dim*time_horizon);
    % Corresponding xmax
    xmax = opt_input_vec_and_xmax(sys.input_dim * time_horizon + 1: end);
end
