function [underapprox_stoch_reach_avoid_polytope, ...
          opt_input_vector_at_vertices, ...
          varargout] = ...
          getUnderapproxStochReachAvoidSet(sys, ...
                                           target_tube, ...
                                           init_safe_set_affine_const, ...
                                           prob_thresh_of_interest, ...
                                           set_of_direction_vectors,...
                                           method,...
                                           varargin)
% SReachTools/stochasticReachAvoid/getUnderapproxStochReachAvoidSet: Obtain an
% open-loop controller-based underaproximative stochastic reach-avoid set using
% Fourier transform, convex optimization, and patternsearch
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
%                                target_tube, ...
%                                init_safe_set_affine_const, ...
%                                prob_thresh_of_interest, ...
%                                set_of_direction_vectors,...
%                                varargin)
% 
% Inputs:
% -------
%   sys                  - LtiSystem object describing the system to be verified
%   target_tube          - Target tube to stay within [TargetTube object]
%   init_safe_set_affine_const        
%                        - Affine constraints (if any) on the initial state
%                          Must include a translate of the affine hull of the
%                          set_of_direction_vectors                          
%   prob_thresh_of_interest 
%                        - Probability thresh (\theta) that defines the
%                          stochastic reach-avoid set 
%                          {x_0: V_0^\ast( x_0) \geq \theta}
%   set_of_direction_vectors
%                        - Number of unique directions defining the polytope
%                          vertices. Its span is the affine hull whose slice of
%                          the stochastic reach-avoid set is of interest.
%   method               - TODO
%   desired_accuracy     - (Optional for 'genzps') Accuracy expected for the
%                          integral of the Gaussian random vector X over the
%                          concatenated_target_tube [Default 5e-3]
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
%                          set_of_direction_vectors
%   R                    - (Optional for ccc) Chebyshev radius associated with
%                          xmax
%   artificial_consv_pwl - (Optional for ccc) Artificial conservativeness due to
%                          the use of piecewise linear approximation
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

    % DEPENDENCY CHECK: Check if dependencies have been installed correctly
    assert(exist('mpt_init','file')==2, ...
           'SReachTools:setup_error', ...
           ['This function uses MPT3. Please get it from ', ...
            'http://control.ee.ethz.ch/~mpt/3/.']);
    assert(exist('cvx_begin','file')==2, ...
           'SReachTools:setup_error', ...
           'This function uses CVX. Please get it from http://cvxr.com.');
    assert(exist('normcdf','file')==2, ...
           'SReachTools:setup_error', ...
           ['This function needs MATLAB''s Statistics and Machine Learning', ...
            ' Toolbox.']);

    
    % Target tubes has polyhedra T_0, T_1, ..., T_{time_horizon}
    time_horizon = length(target_tube)-1;

    % Get half space representation of the target tube and time horizon
    [concat_target_tube_A, concat_target_tube_b] = target_tube.concat();
    n_ineq_init_set = size(target_tube(1).H,1);
    concat_target_tube_A_rand =...
        concat_target_tube_A(n_ineq_init_set+1:end, sys.state_dim+1:end);
    concat_target_tube_b_rand =...
        concat_target_tube_b(n_ineq_init_set+1:end);
    % Compute M --- the number of polytopic halfspaces to worry about
    n_lin_const = size(concat_target_tube_A_rand,1);

    % Construct the constrained initial safe set
    init_safe_set = Polyhedron('H', target_tube(1).H,...
        'He',[target_tube(1).He;init_safe_set_affine_const]);
    
    % Construct U^N 
    % GUARANTEES: Non-empty input sets (polyhedron) and scalar
    %             time_horizon>0
    [concat_input_space_A, concat_input_space_b] = ...
        getConcatInputSpace(sys, time_horizon);

    % Check prob_thresh_of_interest is a scalar in [0.5,1]
    assert( isscalar(prob_thresh_of_interest) &&...
                prob_thresh_of_interest >= 0.5 &&...
                prob_thresh_of_interest <= 1, ...
            'SReachTools:invalidArgs', ...
            'prob_thresh_of_interest must be a scalar in [0.5,1]');

    % Get the n_direction vectors
    assert( size(set_of_direction_vectors, 1) == sys.state_dim, ...
            'SReachTools:invalidArgs', ...
            ['set_of_direction_vectors should be a collection of ',...
             ' n-dimensional column vectors.']);
    assert( size(set_of_direction_vectors, 2) >= 2, ...
            'SReachTools:invalidArgs', ...
            'set_of_direction_vectors should at least have two directions.');
    % TODO: Do input handling for uniqueness
    n_direction_vectors = size(set_of_direction_vectors, 2);

    % Compute H, mean_X_sans_input_sans_initial_state, cov_X_sans_input, and Z.
    % See @LtiSystem\getConcatMats for the notation.
    % GUARANTEES: Gaussian-perturbed LTI system (sys)
    [H, mean_X_sans_input_sans_initial_state, cov_X_sans_input, Z] = ...
        getHmatMeanCovForXSansInput(sys,zeros(sys.state_dim,1),time_horizon);

    switch(lower(method))
        case 'ccc'
            %% Computation of xmax and the associated opt open-loop controller
            % TODO: Translate desired_accuracy to piecewise_count
            [invcdf_approx_m, invcdf_approx_c, lb_deltai, max_err_pwl] =...
                computeNormCdfInvOverApprox();
            % Using PWL approximation introduces an artificial conservativeness.
            % This artificial conservativeness is proportional to the number of
            % linear constraints and the maximum of the the overapproximation
            % error (max_err_pwl) or the closest we can get to zero for
            % normcdfinv (lb_deltai). 
            artificial_consv_pwl = max(max_err_pwl,lb_deltai) * n_lin_const;

            %% PWL code preps
            % Compute \sqrt{h_i^\top * \Sigma_X_no_input * h_i}
            sigma_vector = diag(sqrt(diag(concat_target_tube_A_rand *...
                cov_X_sans_input * concat_target_tube_A_rand')));
            
            %% Solve the optimization problem to compute a point that is deepest
            %% inside the desired set --- First find the maximum value, then
            %% search among it the deepest
            % Maximum value for the open-loop-based control
            % maximize W_0(x,U) 
            % subject to
            %   W_0(x,U) >= prob_thresh_of_interest 
            %                      [Computation of W_0(x,U) is done via
            %                       chance-constraint reformulation and risk
            %                       allocation]
            %   x\in InitState
            %   U\in U^N
            cvx_clear
            cvx_begin quiet
                variable mean_X(sys.state_dim * time_horizon, 1);
                variable deltai(n_lin_const, 1);
                variable norminvover(n_lin_const, 1);
                variable xmax(sys.state_dim, 1);
                variable Umax(sys.input_dim * time_horizon,1);
            
                maximize (1-sum(deltai))
                subject to
                    % Chance-constraint reformulation of the safety prob (1-6)
                    % (1) Slack variables that are bounded below by the
                    % piecewise linear approximation of \phi^{-1}(1-\delta)
                    for deltai_indx=1:n_lin_const
                        norminvover(deltai_indx) >= invcdf_approx_m.*...
                            deltai(deltai_indx) + invcdf_approx_c; 
                    end
                    % (2) Trajectory
                    mean_X == Z * xmax + H * Umax +...
                        mean_X_sans_input_sans_initial_state;
                    % (3) Ono's type of reformulation of chance constraints
                    concat_target_tube_A_rand* mean_X +...
                        sigma_vector * norminvover...  
                            <= concat_target_tube_b_rand;
                    % (4) Lower bound on delta due to pwl's domain
                    deltai >= lb_deltai;
                    % (5) Upper bound on delta to ensure \phi^{-1}(1-\delta) is
                    % concave
                    deltai <= 0.5;
                    % (6) W_0(x,U) >= alpha
                    1-sum(deltai) >= prob_thresh_of_interest
                    % Safety constraints
                    init_safe_set.A * xmax  <= init_safe_set.b
                    init_safe_set.Ae * xmax == init_safe_set.be
                    % Input constraints
                    concat_input_space_A * Umax <= concat_input_space_b;
            cvx_end
            if strcmpi(cvx_status, 'Solved')
                max_underapprox_reach_avoid_prob = 1-sum(deltai);
                %% Compute the largest radius possible. This will ensure that
                %% W_0(x,U) can be given more importance than the radius in the next
                %% optimization problem
                % maximize R 
                % subject to
                %   W_0(x,U) >= W_0(x_max,U_max)
                %   x\in InitState and and admits a radius R ball fit inside the
                %       inequalities
                %   U\in U^N
                % Dual norm of a_i in init_safe_set for chebyshev-centering
                % constraint
                dual_norm_of_init_safe_set_A =...
                    sqrt(diag(init_safe_set.A*init_safe_set.A')); 
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
                        for deltai_indx=1:n_lin_const
                            norminvover(deltai_indx) >= invcdf_approx_m.*...
                                deltai(deltai_indx) + invcdf_approx_c; 
                        end
                        % (2) Trajectory
                        mean_X == Z * xmax + H * Umax +...
                            mean_X_sans_input_sans_initial_state;
                        % (3) Ono's type of reformulation of chance constraints
                        concat_target_tube_A_rand* mean_X +...
                            sigma_vector * norminvover...
                                <= concat_target_tube_b_rand;
                        % (4) Lower bound on delta due to pwl's domain
                        deltai >= lb_deltai;
                        % (5) Upper bound on delta to ensure \phi^{-1}(1-\delta) is
                        % concave
                        deltai <= 0.5;
                        % (6) W_0(x,U) >= W_0(x_max,U_max)
                        1-sum(deltai) >= max_underapprox_reach_avoid_prob
                        % Chebyshev-centering constraints-based safety
                        init_safe_set.A * xmax + ...
                                R * dual_norm_of_init_safe_set_A...
                                    <= init_safe_set.b
                        init_safe_set.Ae * xmax == init_safe_set.be
                        % Input constraints
                        concat_input_space_A * Umax <= concat_input_space_b;
                cvx_end
                opt_input_vector_for_xmax = Umax;
                % Now that we have computed xmax, make sure that the
                % user-provided set_of_direction_vectors will span the affine
                % hull.  Compute xmax + directions \forall directions
                check_set_of_direction_vectors = ...
                    repmat(xmax,1,n_direction_vectors)...
                        + set_of_direction_vectors;
                % Check if all of these vectors satisfy the equality constraints
                % of init_safe_set
                xmaxPlusdirs_within_init_safe_set_eq =...
                    abs(init_safe_set.Ae * check_set_of_direction_vectors -...
                        init_safe_set.be)<eps;
                assert(all(all(xmaxPlusdirs_within_init_safe_set_eq)),...
                    'SReachTools:invalidArgs', ...
                    ['set_of_direction_vectors+xmax should lie in the',...
                     '  affine hull intersecting safe_set.']);
            else
                % Failed to find a single point within the polytope (if there is
                % at least one x for which W(x,U)= prob_thresh_of_interest, then
                % R would have been >=0, and the above problem would have been
                % feasible)
                % % Option 1: Do nothing and return NaN
                % max_underapprox_reach_avoid_prob = NaN;
                % opt_input_vector_for_xmax = nan(sys.input_dim * time_horizon,1);    
                % Option 2: Compute W_max
                disp('Computing W_max');
                cvx_begin quiet
                    variable mean_X(sys.state_dim * time_horizon, 1);
                    variable deltai(n_lin_const, 1);
                    variable norminvover(n_lin_const, 1);
                    variable xmax(sys.state_dim, 1);
                    variable Umax(sys.input_dim * time_horizon,1);

                    maximize (1-sum(deltai))
                    subject to
                        % Chance-constraint reformulation of the safety prob (1-6)
                        % (1) Slack variables that are bounded below by the
                        % piecewise linear approximation of \phi^{-1}(1-\delta)
                        for deltai_indx=1:n_lin_const
                            norminvover(deltai_indx) >= invcdf_approx_m.*...
                                deltai(deltai_indx) + invcdf_approx_c; 
                        end
                        % (2) Trajectory
                        mean_X == Z * xmax + H * Umax +...
                            mean_X_sans_input_sans_initial_state;
                        % (3) Ono's type of reformulation of chance constraints
                        concat_target_tube_A_rand* mean_X +...
                            sigma_vector * norminvover...  
                                <= concat_target_tube_b_rand;
                        % (4) Lower bound on delta due to pwl's domain
                        deltai >= lb_deltai;
                        % (5) Upper bound on delta to ensure \phi^{-1}(1-\delta) is
                        % concave
                        deltai <= 0.5;
                        % We can't go below 0.5 with this method
                        (1-sum(deltai)) >= 0.5;
                        % Safety constraints
                        init_safe_set.A * xmax  <= init_safe_set.b
                        init_safe_set.Ae * xmax == init_safe_set.be
                        % Input constraints
                        concat_input_space_A * Umax <= concat_input_space_b;
                cvx_end
                if strcmpi(cvx_status, 'Solved')
                    max_underapprox_reach_avoid_prob = 1-sum(deltai);
                    opt_input_vector_for_xmax = Umax;    
                else
                    max_underapprox_reach_avoid_prob = NaN;
                    opt_input_vector_for_xmax = nan(sys.input_dim * time_horizon,1);                    
                end
            end
            % If non-trivial underapproximative stochastic reach-avoid polytope
            if max_underapprox_reach_avoid_prob < prob_thresh_of_interest || isnan(max_underapprox_reach_avoid_prob)
                % Stochastic reach-avoid underapproximation is empty and no
                % admissible open-loop policy exists
                fprintf(['Polytopic underapproximation does not exist for ',...
                         'alpha = %1.2f since maximum reach probability ',...
                         '(%1.3f) with open-loop is below alpha.\n\n'], ...
                         prob_thresh_of_interest,...
                         max_underapprox_reach_avoid_prob);
                % Assigning the outputs to trivial results
                underapprox_stoch_reach_avoid_polytope = Polyhedron();
                xmax = nan(sys.state_dim, 1);
                opt_input_vector_at_vertices =...
                    nan(sys.input_dim * time_horizon, n_direction_vectors);
                opt_theta_i = nan(1, n_direction_vectors);
                opt_reachAvoid_i = nan(1, n_direction_vectors);
                vertex_underapprox_polytope = nan(sys.state_dim,...
                    n_direction_vectors);
                R = nan;
            else
                % Stochastic reach-avoid underapproximation is non-trivial
%                 fprintf(['Polytopic underapproximation exists for alpha = ',...
%                          '%1.2f since W(x_max) = %1.3f.\n\n'], ...
%                          prob_thresh_of_interest, ...
%                          max_underapprox_reach_avoid_prob);

                % Pre-allocation of relevant outputs
                opt_input_vector_at_vertices =...
                    zeros(sys.input_dim * time_horizon,n_direction_vectors);
                opt_theta_i = zeros(1, n_direction_vectors);
                opt_reachAvoid_i = zeros(1, n_direction_vectors);

                %% Iterate over all direction vectors + xmax
                for direction_index = 1: n_direction_vectors
                    % Get direction_index-th direction in the hyperplane
                    direction = set_of_direction_vectors(:,direction_index);

                    fprintf('Analyzing direction :%2d/%2d\n',...
                        direction_index, n_direction_vectors);

                    %% Solve the optimization problem to compute the boundary
                    %% point of the desired underapproximative stochastic
                    %% reach-avoid set
                    % maximize theta
                    % subject to
                    %   boundary_point = xmax + theta * direction
                    %   U\in U^N
                    %   boundary_point\in InitState
                    %   W_0( boundary_point, U) >= prob_thresh_of_interest 
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
                            boundary_point == xmax + theta *direction
                            % Safe boundary point
                            init_safe_set.A * boundary_point <= init_safe_set.b
                            init_safe_set.Ae* boundary_point == init_safe_set.be
                            % Input constraints
                            concat_input_space_A*U_vector<=concat_input_space_b;
                            % Chance-constraint reformulation of safety
                            % prob(1-6) (1) Slack variables that are bounded
                            % below by the piecewise linear approximation of
                            % \phi^{-1}(1-\delta)
                            for deltai_indx=1:n_lin_const
                                norminvover(deltai_indx) >= invcdf_approx_m.*...
                                    deltai(deltai_indx) + invcdf_approx_c; 
                            end
                            % (2) Trajectory
                            mean_X == Z * boundary_point + H * U_vector +...
                                mean_X_sans_input_sans_initial_state;
                            % (3) Ono's type of reformulation of chance
                            % constraints
                            concat_target_tube_A_rand * mean_X +...
                                sigma_vector * norminvover...
                                    <= concat_target_tube_b_rand;
                            % (4) Lower bound on delta due to pwl's domain
                            deltai >= lb_deltai;
                            % (5) Upper bound on delta to ensure
                            % \phi^{-1}(1-\delta) is concave
                            deltai <= 0.5;
                            % (6) W_0(x,U) >= alpha
                            1-sum(deltai) >= prob_thresh_of_interest
                    cvx_end
                    opt_theta_i(direction_index) = theta;
                    opt_input_vector_at_vertices(:,direction_index) = U_vector;
                    opt_reachAvoid_i(direction_index) = 1 - sum(deltai);
                    vertex_underapprox_polytope(:, direction_index)=...
                        boundary_point;
                end
                underapprox_stoch_reach_avoid_polytope = Polyhedron('V', ...
                    vertex_underapprox_polytope');
            end
        case 'genzps'
            % Parsing the optional arguments 
            if length(varargin) == 2
                % First optional argument is the desired_accuracy
                assert(isscalar(varargin{1}), ...
                       'SReachTools:invalidArgs', ...
                       'Expected a scalar value for desired_accuracy');
                desired_accuracy = varargin{1};
                % Second optional argument is the options for patternsearch,
                % PSoptions (TODO: No validation being done here)
                PSoptions = varargin{2};
            elseif isempty(varargin)
                display_string = 'off';        % Turned it off to use the
                                               % aligned output
                desired_accuracy = 1e-3;
                PSoptions = psoptimset('Display', display_string);
            else
                exc = SrtInvalidArgsError(['desired_accuracy and ', ...
                    'PSoptions together are the only additional options.']);
                throwAsCaller(exc);
            end

            %% Computation of xmax and the associated opt open-loop controller
            disp(['Computing the x_max for the Fourier transform-based ', ...
                  'underapproximation']);
            [max_underapprox_reach_avoid_prob, ...
             xmax, ...
             opt_input_vector_for_xmax] = ...
                computeXmaxForStochReachAvoidSetUnderapprox(...
                    sys, ...
                    time_horizon, ...
                    init_safe_set, ...
                    concat_input_space_A, ... 
                    concat_input_space_b, ...
                    concat_target_tube_A_rand, ... 
                    concat_target_tube_b_rand, ...
                    Z, ...
                    H, ...
                    mean_X_sans_input_sans_initial_state, ...
                    cov_X_sans_input, ...
                    desired_accuracy, ...
                    PSoptions);

            % If non-trivial underapproximative stochastic reach-avoid polytope
            if max_underapprox_reach_avoid_prob <...
                prob_thresh_of_interest
                % Stochastic reach-avoid underapproximation is empty and no
                % admissible open-loop policy exists
                fprintf(['Polytopic underapproximation does not exist for ',...
                         'alpha = %1.2f since W(x_max) = %1.3f.\n\n'], ...
                         prob_thresh_of_interest, ...
                         max_underapprox_reach_avoid_prob);

                % Assigning the outputs to trivial results
                underapprox_stoch_reach_avoid_polytope = Polyhedron();
                opt_input_vector_at_vertices = nan(...
                    sys.input_dim * time_horizon, ...
                    n_direction_vectors);
                opt_theta_i = zeros(1, n_direction_vectors);
                opt_reachAvoid_i = zeros(1, n_direction_vectors);
                vertex_underapprox_polytope = nan(sys.state_dim,...
                    n_direction_vectors);
            else
                % Stochastic reach-avoid underapproximation is non-trivial
                fprintf(['Polytopic underapproximation exists for alpha = ',...
                         '%1.2f since W(x_max) = %1.3f.\n\n'], ...
                         prob_thresh_of_interest, ...
                         max_underapprox_reach_avoid_prob);

                % For storing boundary points
                opt_theta_i = zeros(1, n_direction_vectors);
                opt_reachAvoid_i = zeros(1, n_direction_vectors);
                opt_input_vector_at_vertices = nan(...
                    sys.input_dim * time_horizon, ...
                    n_direction_vectors);

                %% Iterate over all direction vectors + xmax
                for direction_index = 1: n_direction_vectors
                    % Get direction_index-th direction in the hyperplane
                    direction = set_of_direction_vectors(:,direction_index);

                    fprintf('Analyzing direction :%2d/%2d\n',...
                        direction_index, n_direction_vectors);

                    %% Bounds on theta \in [lb_on_theta, ub_on_theta] 
                    %% Lower bound is always 0 as xmax could be a vertex
                    lb_theta = 0;
                    tol_bisection = 1e-2;  %TODO
                    % Computation of the upper bound
                    A_times_direction = init_safe_set.A * direction;
                    %% Compute theta_max for the given direction and update
                    %% ub_theta for each direction, given by the optimization
                    %% problem
                    % minimize -theta
                    % s.t.    theta*(A_init_safe_set*direction) <=
                    %               init_safe_set.b-init_safe_set.A*xmax
                    cvx_begin quiet
                        variable theta(1)
                        minimize -theta
                        subject to
                            theta*A_times_direction <= init_safe_set.b...
                                                           -init_safe_set.A*xmax
                    cvx_end
                    ub_theta = theta;
                    fprintf('\b | Upper bound of theta: %1.2f\n',ub_theta);

                    %% Bisection-based computation of the maximum extension of
                    %%the ray originating from xmax
                    fprintf(['OptRAProb | OptTheta | LB_theta | UB_theta | ',...
                                 'OptInp^2 | Exit reason\n']); %10 characters b/n | |
                    
                    % Temp variables that are updated as part of bisection
                    opt_inputs_so_far = opt_input_vector_for_xmax;
                    opt_reachAvoid_prob_so_far = max_underapprox_reach_avoid_prob;
                    opt_theta_so_far = lb_theta;
                    % While the gap remains above tolerance_bisection
                    while (ub_theta - lb_theta) > tol_bisection
                        % Set phi as the middle point of the interval
                        phi = (ub_theta+lb_theta)/2;
                        % Candidate initial state provides mean_X_sans_input
                        % See @LtiSystem/getConcatMats for notation
                        candidate_initial_state = xmax + phi * direction;
                        mean_X_sans_input = Z * candidate_initial_state + ...
                            mean_X_sans_input_sans_initial_state;
                        % Compute the maximum reach-avoid prob and the
                        % associated open-loop opt controller starting from
                        % candidate_initial_state using Genz+patternsearch
                        % Add desired_accuracy to prob_thresh_of_interest so that we
                        % are guaranteed to meet the desired thresh 
                        negLogReachAvoidProbGivenInputVectorFeas= @(input_vector)...
                            -log(min(computeReachAvoidProb(input_vector, ...
                                        mean_X_sans_input, ...
                                        cov_X_sans_input, ...
                                        H, ...
                                        concat_target_tube_A_rand, ...
                                        concat_target_tube_b_rand, ...
                                        desired_accuracy),...
                                     prob_thresh_of_interest+5*desired_accuracy));

                        % Compute the opt admissible input_vector that ensures
                        % W_0(x,U)\geq \alpha, given x_0
                        %
                        % minimize
                        % -log(min(ReachAvoidProb(U),prob_thresh_of_interest))
                        % subject to
                        %        concat_input_space_A * U \leq concat_input_space_b
                        [opt_input_at_this_step, patternsearch_opt_value]= ...
                              patternsearch(...
                                negLogReachAvoidProbGivenInputVectorFeas,...
                                opt_inputs_so_far, ...
                                concat_input_space_A, ...
                                concat_input_space_b, ...
                                [],[],[],[],[], ...
                                PSoptions);
                        prob_value_at_this_step = exp(-patternsearch_opt_value);
                        if  prob_value_at_this_step -...
                                prob_thresh_of_interest >= desired_accuracy/2
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
                        % Print current solution
                        fprintf(strcat(...
                            '  %1.4f  |  %1.4f  |  %1.4f  |  %1.4f  |  %1.4f  |  ', ...
                            exitmessage, ' \n'), ...
                                   opt_reachAvoid_prob_so_far, ...
                                   opt_theta_so_far, ...
                                   lb_theta, ...
                                   ub_theta, ...
                                   opt_inputs_so_far'*opt_inputs_so_far);

                        % Update bounds
                        if  strcmpi(exitmessage, ' Feasible')
                            % Theta proved to be feasible
                            lb_theta = phi;
                        else
                            % Shrink the upper bound
                            ub_theta = phi;
                        end
                    end        
                    opt_theta_i(direction_index) = opt_theta_so_far;
                    opt_input_vector_at_vertices(:,direction_index) = opt_inputs_so_far;
                    opt_reachAvoid_i(direction_index) = opt_reachAvoid_prob_so_far;
                end
                vertex_underapprox_polytope = xmax +...
                                      opt_theta_i.*set_of_direction_vectors;
                %% Construction of underapprox_stoch_reach_avoid_polytope
                underapprox_stoch_reach_avoid_polytope = Polyhedron('V', ...
                    vertex_underapprox_polytope');
            end
        otherwise
            exc = SrtInvalidArgsError(['Unsupported method requested in ', ...
                   'getUnderapproxStochReachAvoidSet.']);
            throwAsCaller(exc);
    end
    varargout{1} = xmax;
    varargout{2} = opt_input_vector_for_xmax;
    varargout{3} = max_underapprox_reach_avoid_prob;
    varargout{4} = opt_theta_i;
    varargout{5} = opt_reachAvoid_i;
    varargout{6} = vertex_underapprox_polytope;
    if strcmpi(method, 'ccc')
        varargout{7} = R;
        varargout{8} = artificial_consv_pwl;
    end
end


%method = 'best'; % 'fast'/'best'
%tol_bisection = 1e-2;
%desired_accuracy = 1e-2;
%PSoptions = psoptimset('Display','off');
%
%%% Compute theta_upper for the given direction
%cvx_begin quiet
    %variable theta_upper nonnegative

    %maximize theta_upper
    %subject to
        %init_safe_set.A * (xmax + theta_upper * direction) <=...
          %init_safe_set.b
        %theta_upper >= theta
%cvx_end
%theta_lower = theta;
%%% Did ccc reach the boundary of init_safe_set or is the method set
%%% to fast (skip genz+ps-based improvement)?
%if (theta_upper - theta_lower) < tol_bisection ||...
    %strcmpi(method,'fast')
    %%% No point doing Genz+ps because ccc found the boundary of
    %%% init_safe_set
    %vertex_underapprox_polytope(:, direction_index)= boundary_point;
%else
    %%% Do bisection in the remaining gap using Genz+ps
    %disp('Chance constraint formulation did not reach boundary');
    
%end
%fprintf(['Analyzing direction (shown transposed) ', ...
%         ':%d/%d\n'],...
%         direction_index,n_direction_vectors);
%disp(direction');
