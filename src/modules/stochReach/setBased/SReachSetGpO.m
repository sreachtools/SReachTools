function varargout = SReachSetGpO(method_str, sys, prob_thresh, safety_tube,...
    options)
% Obtain an open-loop controller-based underaproximative stochastic reach-avoid
% set using Fourier transform, convex optimization, and patternsearch
% =============================================================================
%
% getUnderapproxStochReachAvoidSet computes the open-loop controller-based
% underapproximative stochastic reach-avoid set  to the terminal hitting-time
% stochastic reach-avoid problem discussed in
%
% A. Vinod, and M. Oishi, "Scalable Underapproximative Verification of
% Stochastic LTI Systems Using Convexity and Compactness," in Proceedings of
% Hybrid Systems: Computation and Control (HSCC), 2018. 
%
% Specifically, we use chance-constrained approach to speed up the computation.
%
% USAGE: TODO
%
% =============================================================================
%
% [underapprox_stoch_reach_avoid_polytope, ...
%  opt_input_vector_at_vertices, ...
%  varargout] = ...
%  getUnderapproxStochReachAvoidSet(...
%                                sys, ...
%                                safety_tube, ...
%                                init_safe_set_affine_const, ...
%                                prob_thresh, ...
%                                set_of_dir_vecs,...
%                                varargin)
% 
% Inputs:
% -------
%   sys                  - LtiSystem object describing the system to be verified
%   safety_tube          - Target tube to stay within [Tube object]
%   init_safe_set_affine_const        
%                        - Affine constraints (if any) on the initial state
%                          Must include a translate of the affine hull of the
%                          set_of_dir_vecs                          
%   prob_thresh 
%                        - Probability thresh (\theta) that defines the
%                          stochastic reach-avoid set 
%                          {x_0: V_0^\ast( x_0) \geq \theta}
%   set_of_dir_vecs
%                        - Number of unique directions defining the polytope
%                          vertices. Its span is the affine hull whose slice of
%                          the stochastic reach-avoid set is of interest.
%   method               - TODO
%   desired_accuracy     - (Optional for 'genzps') Accuracy expected for the
%                          integral of the Gaussian random vector X over the
%                          concatenated_safety_tube [Default 5e-3]
%   PSoptions            - (Optional for 'genzps') Options for patternsearch 
%                          [Default psoptimset('Display', 'off')]
%
% Outputs:
% --------
%   underapprox_stoch_reach_avoid_polytope
%                        - Underapproximative polytope of dimension
%                          sys.state_dim which underapproximates the
%                          terminal-hitting stochastic reach avoid set
%   opt_input_vector_at_vertices 
%                        - Optimal open-loop policy ((sys.input_dim) *
%                          time_horizon)-dim.  vector U = [u_0; u_1; ...; u_N]
%                          (column vector) for each vertex of the polytope
%   xmax                 - (Optional) Initial state that has the maximum
%                          stochastic reach-avoid prob using an open-loop
%                          controller
%   opt_input_vector_for_xmax
%                        - (Optional) Optimal open-loop policy
%                          ((sys.input_dim) * time_horizon)-dimensional
%                          vector U = [u_0; u_1; ...; u_N] (column vector) for
%                          xmax
%   max_underapprox_reach_avoid_prob
%                        - (Optional) Maximum attainable stochastic reach-avoid
%                          prob using an open-loop controller; Maximum
%                          terminal-hitting time reach-avoid prob at xmax
%   opt_theta_i          - (Optional) Vector comprising of scaling factors along
%                          each direction of interest
%   opt_reachAvoid_i     - (Optional) Maximum terminal-hitting time reach-avoid
%                          prob at the vertices of the polytope
%   vertex_underapprox_polytope
%                        - (Optional) Vertices of the polytope: xmax + opt_theta_i *
%                          set_of_dir_vecs
%   R                    - (Optional for ccc) Chebyshev radius associated with
%                          xmax
%
% See also examples/FtCVXUnderapproxVerifyCWH.mlx*.
%
% Notes:
% ------
% * NOT ACTIVELY TESTED: Builds on other tested functions.
% * MATLAB DEPENDENCY: Uses MATLAB's Global Optimization Toolbox; Statistics and
%                      Machine Learning Toolbox.
%                      Needs patternsearch for gradient-free optimization
%                      Needs normpdf, normcdf, norminv for Genz's algorithm
% * EXTERNAL DEPENDENCY: Uses MPT3 and CVX
%                      Needs MPT3 for defining a controlled system and the
%                      definition of the safe, the target (polytopic) sets, and
%                      the affine hull of interest
%                      Needs CVX to setup convex optimization problems that
%                      1) initializes the patternsearch-based optimization, and
%                      2) computes the upper bound for the bisection
% * Specify both desired_accuracy and PSoptions or neither to use the defaults 
% * max_underapprox_reach_avoid_prob is the highest thresh
%   that may be given while obtaining a non-trivial underapproximation
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

    validatestring(method_str,{'genzps-open'});
    
    inpar = inputParser();
    inpar.addRequired('sys', @(x) validateattributes(x,...
        {'LtiSystem','LtvSystem'}, {'nonempty'}));
    inpar.addRequired('prob_thresh', @(x) validateattributes(x, {'numeric'},...
        {'scalar','>=',0,'<=',1}));
    inpar.addRequired('safety_tube',@(x) validateattributes(x,{'Tube'},...
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
    [concat_input_space_A, concat_input_space_b] = getConcatInputSpace(sys,...
        time_horizon);
    
    % Compute the concatenated state dynamic matrices
    [Z, H, ~] = getConcatMats(sys, time_horizon); 
    
    % Compute the mean of G*W vector --- mean of concatenated state vector
    % under zero input and zero state conditions
    sysnoi = LtvSystem('StateMatrix',sys.state_mat,'DisturbanceMatrix',...
        sys.dist_mat,'Disturbance',sys.dist);
    [mean_X_zizs, cov_X_sans_input] = SReachFwd('concat-stoch', sysnoi,...
        zeros(sys.state_dim,1), time_horizon);

    % Construct the constrained initial safe set
    init_safe_set = Polyhedron('H', safety_tube(1).H,...
                      'He',[safety_tube(1).He;options.init_safe_set_affine.He]);
    
    %% Step 1: Find xmax
    % Compute the chebyshev center of the initial safe set and seek the
    % optimal xmax that maximizes the open-loop safety from there
    dual_norm_of_init_safe_set_A = sqrt(diag(init_safe_set.A*init_safe_set.A')); 
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
    initial_guess_input_vector_and_xmax = [guess_concatentated_input_vector;
                                           initial_x_for_xmax];
    
    % Construct the reach-avoid cost function: -log(ReachAvoidProbability(x_0,U))
    % input_vector_and_xmax has first the unrolled input vector followed by
    % the initial state that is optimized
    negLogReachProbGivenInputVecInitState = @(input_vector_and_xmax)...
      -log(computeReachAvoidProb(...
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

    % Compute xmax, the input policy, and the max reach-avoid probability
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
    % Maximum attainable terminal hitting-time stochastic reach-avoid
    % probability using open-loop controller
    max_underapprox_reach_prob = exp(-opt_neg_log_reach_prob);
    % Optimal open_loop_control_policy
    opt_input_vec_for_xmax=opt_input_vec_and_xmax(1:sys.input_dim*time_horizon);
    % Corresponding xmax
    xmax = opt_input_vec_and_xmax(sys.input_dim * time_horizon + 1: end);

    % No. of directions
    n_dir_vecs = size(options.set_of_dir_vecs, 2);
    
    % If non-trivial underapproximative stochastic reach-avoid polytope
    if max_underapprox_reach_prob < prob_thresh
        % Stochastic reach-avoid underapproximation is empty and no
        % admissible open-loop policy exists
        if options.verbose >= 1
            fprintf(['Polytopic underapproximation does not exist for ',...
                 'alpha = %1.2f since W(x_max) = %1.3f.\n\n'], ...
                 prob_thresh, max_underapprox_reach_prob);
        end

        % Assigning the outputs to trivial results
        underapprox_stoch_reach_polytope = Polyhedron();
        opt_input_vec_at_vertices =nan(sys.input_dim * time_horizon,n_dir_vecs);
        opt_theta_i = zeros(1, n_dir_vecs);
        opt_reach_prob_i = zeros(1, n_dir_vecs);
        vertex_underapprox_polytope = nan(sys.state_dim, n_dir_vecs);
    else
        if options.verbose >= 1
            % Stochastic reach-avoid underapproximation is non-trivial
            fprintf(['Polytopic underapproximation exists for alpha = ',...
                     '%1.2f since W(x_max) = %1.3f.\n\n'], ...
                     prob_thresh, max_underapprox_reach_prob);
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
            ub_theta = theta;
            if options.verbose >= 1
                fprintf('\b | Upper bound of theta: %1.2f\n',ub_theta);

                fprintf(['OptRAProb |  OptTheta  |  LB_theta  |  UB_theta  ',...
                          '|  OptInp^2  | Exit reason\n']); %10 characters b/n | |
            end
            %% Bisection-based computation of the maximum extension of
            %% the ray originating from xmax                
            % Temp variables that are updated as part of bisection
            opt_inputs_so_far = opt_input_vec_for_xmax;
            opt_reachAvoid_prob_so_far = max_underapprox_reach_prob;
            opt_theta_so_far = lb_theta;
            % While the gap remains above tolerance_bisection
            while (ub_theta - lb_theta) > options.tol_bisect
                % Set phi as the middle point of the interval
                phi = (ub_theta+lb_theta)/2;
                % Candidate initial state provides mean_X_sans_input
                % See @LtiSystem/getConcatMats for notation
                candidate_initial_state = xmax + phi * dir_vec;
                mean_X_sans_input = Z * candidate_initial_state + mean_X_zizs;
                % Compute the maximum reach-avoid prob and the
                % associated open-loop opt controller starting from
                % candidate_initial_state using Genz+patternsearch
                % Add desired_accuracy to prob_thresh so that we
                % are guaranteed to meet the desired thresh 
                negLogReachProbGivenInputVectorFeas= @(input_vector)...
                    -log(min(computeReachAvoidProb(input_vector, ...
                                mean_X_sans_input, ...
                                cov_X_sans_input, ...
                                H, ...
                                concat_safety_tube_A, ...
                                concat_safety_tube_b, ...
                                options.desired_accuracy),...
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
                        negLogReachProbGivenInputVectorFeas,...
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
                    exitmessage = sprintf(' Infeasible (%1.3f)',...
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
        vertex_underapprox_polytope = xmax+opt_theta_i.*options.set_of_dir_vecs;
        %% Construction of underapprox_stoch_reach_avoid_polytope
        underapprox_stoch_reach_polytope = Polyhedron('V', ...
            vertex_underapprox_polytope');
    end
    varargout{1} = underapprox_stoch_reach_polytope;
    extra_info.xmax = xmax;
    extra_info.opt_input_vec_for_xmax = opt_input_vec_for_xmax;
    extra_info.opt_input_vec_at_vertices = opt_input_vec_at_vertices;
    extra_info.max_underapprox_reach_prob = max_underapprox_reach_prob;
    extra_info.opt_theta_i = opt_theta_i;
    extra_info.opt_reach_prob_i = opt_reach_prob_i;
    extra_info.vertex_underapprox_polytope = vertex_underapprox_polytope;
    if nargout > 1
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
        throwAsCaller(SrtInvalidArgsError(['Mismatch in prob_str in the ',...
            'options']));
    end
    if ~strcmpi(options.method_str,method_str)
        throwAsCaller(SrtInvalidArgsError(['Mismatch in method_str in the ',...
            'options']));
    end        
    
    % Make sure the user specified set_of_dir_vecs is valid
    if size(options.set_of_dir_vecs, 1) ~= sys.state_dim
        throwAsCaller(SrtInvalidArgsError(['set_of_dir_vecs should be a ',...
            'collection of n-dimensional column vectors.']));
    end
    if size(options.set_of_dir_vecs, 2) < 2, ...
        throwAsCaller(SrtInvalidArgsError(['set_of_dir_vecs should at least',...
            ' have two directions.']));
    end
    % Make sure the user specified init_safe_set_affine is of the correct
    % dimension
    if options.init_safe_set_affine.Dim ~= sys.state_dim &&...
            ~isempty(options.init_safe_set_affine.H)
        throwAsCaller(SrtInvalidArgsError(['init_safe_set_affine must be ',...
            'an sys.state_dim-dimensional affine set']));
    end    
end
