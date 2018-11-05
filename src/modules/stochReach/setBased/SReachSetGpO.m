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
%
% Notes:
% ------
% * extra_info.xmax_reach_prob is the highest prob_thresh that may be given
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
    
    % Half-space representation for the safety tube
    [concat_safety_tube_A, concat_safety_tube_b] =...
        safety_tube.concat([1 time_horizon]+1);
    
    % Half-space representation for the concatenated input space
    [concat_input_space_A, concat_input_space_b] = getConcatInputSpace(sys, ...
        time_horizon);
    
    % Compute the concatenated state dynamic matrices
    [Z, H, ~] = getConcatMats(sys, time_horizon); 
    
    % Compute the mean of G*W vector --- mean of concatenated state vector
    % under zero input and zero state conditions
    sysnoi = LtvSystem('StateMatrix',sys.state_mat,'DisturbanceMatrix', ...
        sys.dist_mat,'Disturbance',sys.dist);
    [mean_X_zizs, cov_X_sans_input] = SReachFwd('concat-stoch', sysnoi, ...
        zeros(sys.state_dim,1), time_horizon);

    % Construct the constrained initial safe set
    init_safe_set = Polyhedron('H', safety_tube(1).H, ...
                      'He',[safety_tube(1).He;options.init_safe_set_affine.He]);
    
    %% Step 1: Find xmax
    % Compute the chebyshev center of the initial safe set and seek the
    % optimal xmax that maximizes the open-loop safety from there
    dual_norm_of_init_safe_set_A = norms(init_safe_set.A,2,2); 
    % maximize R - 0.01 |U| 
    % subject to
    %   X = Abar x_0 + H * U + G_matrix * \mu_W 
    %   U \in \mathcal{U}^N
    %   X \in concatenated_safety_tube
    %   x_0\in AffineHull
    %   init_safe_set.A_i * x_0 + R* || init_safe_set.A_i || \leq b_i 
    %                                      (see Boyd's CVX textbook, pg. 418,
    %                                       Chebyshev centering for a polytope)
    cvx_begin quiet
        variable resulting_X_for_xmax(sys.state_dim * time_horizon) 
        variable guess_concatentated_input_vector(sys.input_dim * time_horizon)
        variable initial_x_for_xmax(sys.state_dim)
        variable R

        maximize R - 0.01 * norm(guess_concatentated_input_vector)

        subject to
            R >= 0
            resulting_X_for_xmax == (Z * initial_x_for_xmax ...
                            + H * guess_concatentated_input_vector...
                            + mean_X_zizs)
            concat_input_space_A * guess_concatentated_input_vector <= ...
                                                      concat_input_space_b 
            concat_safety_tube_A * resulting_X_for_xmax <= ...
                            concat_safety_tube_b
            init_safe_set.Ae * initial_x_for_xmax == ...
                                                   init_safe_set.be
            for i = 1:length(init_safe_set.A)
                init_safe_set.A(i,:) * initial_x_for_xmax...
                     + R * dual_norm_of_init_safe_set_A(i) <= init_safe_set.b(i)
            end
    cvx_end
    if ~strcmpi(cvx_status, 'Solved')
        warning('SReachTools:runtime', ['CVX failed to obtain the Chebyshev',...
            'center, potentially due to numerical issues.']); 
    end
    initial_guess_input_vector_and_xmax = [guess_concatentated_input_vector;
                                           initial_x_for_xmax];
    
    % Construct the reach cost function: 
    %   -log(ReachAvoidProbability(x_0,U))
    % input_vector_and_xmax has first the unrolled input vector followed by
    % the initial state that is optimized
    negLogReachProbGivenInputVecInitState = @(input_vector_and_xmax)...
      -log(computeReachProb(...
                input_vector_and_xmax(1: sys.input_dim * time_horizon), ...
                Z * input_vector_and_xmax(...
                    sys.input_dim * time_horizon + 1: end) + mean_X_zizs, ...
                cov_X_sans_input, ...
                H, ...
                concat_safety_tube_A, ...
                concat_safety_tube_b, ...
                options.desired_accuracy));
    
    % Constraint generation --- decision variable [input_vector;xmax]
    input_vector_augmented_affine_hull_Aeq = ...
        [zeros(size(init_safe_set.Ae,1), ...
                            sys.input_dim * time_horizon), ...
                                                 init_safe_set.Ae];
    input_vector_augmented_affine_hull_beq = init_safe_set.be;
    input_vector_augmented_safe_Aineq = blkdiag(concat_input_space_A, ...
                                                init_safe_set.A);
    input_vector_augmented_safe_bineq = [concat_input_space_b;
                                         init_safe_set.b];

    % Compute xmax, the input policy, and the max reach probability
    [opt_input_vec_and_xmax, opt_neg_log_reach_prob]= ...
              patternsearch(...
                negLogReachProbGivenInputVecInitState, ...
                initial_guess_input_vector_and_xmax, ...
                input_vector_augmented_safe_Aineq, ...
                input_vector_augmented_safe_bineq, ...
                input_vector_augmented_affine_hull_Aeq, ...
                input_vector_augmented_affine_hull_beq, ...
                [],[],[], ...
                options.PSoptions);
    
    %% Parse the output of patternsearch
    % Maximum attainable terminal hitting-time stochastic reach
    % probability using open-loop controller
    xmax_reach_prob = exp(-opt_neg_log_reach_prob);
    % Optimal open_loop_control_policy
    Umax=opt_input_vec_and_xmax(1:sys.input_dim*time_horizon);
    % Corresponding xmax
    xmax = opt_input_vec_and_xmax(sys.input_dim * time_horizon + 1: end);

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
                     '%1.2f since W(x_max) = %1.3f.\n\n'], ...
                     prob_thresh, xmax_reach_prob);
        end

        % For storing boundary points
        opt_theta_i = zeros(1, n_dir_vecs);
        opt_reach_prob_i = zeros(1, n_dir_vecs);
        opt_input_vec_at_vertices =nan(sys.input_dim * time_horizon,n_dir_vecs);

        %% Iterate over all direction vectors + xmax and perform a line search 
        %% via bisection
        for dir_index = 1: n_dir_vecs
            % Get direction_index-th direction in the hyperplane
            dir_vec = options.set_of_dir_vecs(:,dir_index);
            
            if options.verbose >= 1
                fprintf('Analyzing direction :%2d/%2d\n',dir_index, n_dir_vecs);
            end

            % Bounds on theta \in [lb_theta, ub_theta] 
            % Lower bound is always 0 as xmax could be a vertex
            lb_theta = 0;
            
            % Computation of the upper bound
            A_times_dir_vec = init_safe_set.A * dir_vec;
            %% Compute theta_max for the given direction and update
            %% ub_theta for each direction, given by the optimization
            %% problem
            % minimize -theta
            % s.t.   theta*A_times_dir_vec<=init_safe_set.b-init_safe_set.A*xmax
            cvx_begin quiet
                variable theta(1)
                minimize -theta
                subject to
                    theta*A_times_dir_vec<=init_safe_set.b-init_safe_set.A*xmax
            cvx_end
            if ~strcmpi(cvx_status, 'Solved')
                warning('SReachTools:runtime', sprintf(['CVX failed to ',...
                    'obtain the upper bound on theta for direction %d/%d, ',...
                    'potentially due to numerical issues.'], dir_index,...
                    n_dir_vecs); 
            end
            ub_theta = theta;
            if options.verbose >= 1
                fprintf('\b | Upper bound of theta: %1.2f\n',ub_theta);

                fprintf(['OptRAProb |  OptTheta  |  LB_theta  |  UB_theta  ',...
                          '|  OptInp^2  | Exit reason\n']); %10 chars b/n | |
            end
            %% Bisection-based computation of the maximum extension of
            %% the ray originating from xmax                
            % Temp variables that are updated as part of bisection
            opt_inputs_so_far = Umax;
            opt_reachAvoid_prob_so_far = xmax_reach_prob;
            opt_theta_so_far = lb_theta;
            % While the gap remains above tolerance_bisection
            while (ub_theta - lb_theta) > options.tol_bisect
                % Set phi as the middle point of the interval
                phi = (ub_theta+lb_theta)/2;
                % Candidate initial state provides mean_X_sans_input
                % See @LtiSystem/getConcatMats for notation
                candidate_initial_state = xmax + phi * dir_vec;
                mean_X_sans_input = Z * candidate_initial_state + mean_X_zizs;
                % Compute the maximum reach prob and the
                % associated open-loop opt controller starting from
                % candidate_initial_state using Genz+patternsearch
                % Add desired_accuracy to prob_thresh so that we
                % are guaranteed to meet the desired thresh 
                negLogReachProbGivenInputVectorFeas= @(input_vector)...
                    -log(min(computeReachProb(input_vector, ...
                                mean_X_sans_input, ...
                                cov_X_sans_input, ...
                                H, ...
                                concat_safety_tube_A, ...
                                concat_safety_tube_b, ...
                                options.desired_accuracy), ...
                                prob_thresh + 5*options.desired_accuracy));

                % Compute the opt admissible input_vector that ensures
                % W_0(x,U)\geq \alpha, given x_0
                %
                % minimize
                % -log(min(ReachAvoidProb(U),prob_thresh))
                % subject to
                %        concat_input_space_A * U \leq concat_input_space_b
                [opt_input_at_this_step, patternsearch_opt_value]= ...
                      patternsearch(...
                        negLogReachProbGivenInputVectorFeas, ...
                        opt_inputs_so_far, ...
                        concat_input_space_A, ...
                        concat_input_space_b, ...
                        [],[],[],[],[], ...
                        options.PSoptions);
                prob_value_at_this_step = exp(-patternsearch_opt_value);
                if  prob_value_at_this_step - prob_thresh >=...
                        options.desired_accuracy/2
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
        vertices_underapprox_polytope=xmax+opt_theta_i.*options.set_of_dir_vecs;
        %% Construction of underapprox_stoch_reach_avoid_polytope
        polytope = Polyhedron('V', vertices_underapprox_polytope');
    end
    varargout{1} = polytope;
    if nargout > 1
        extra_info.xmax = xmax;
        extra_info.Umax = Umax;
        extra_info.xmax_reach_prob = xmax_reach_prob;
        extra_info.opt_theta_i = opt_theta_i;
        extra_info.opt_input_vec_at_vertices = opt_input_vec_at_vertices;
        extra_info.opt_reach_prob_i = opt_reach_prob_i;
        extra_info.vertices_underapprox_polytope= vertices_underapprox_polytope;
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
