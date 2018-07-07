function [underapprox_stoch_reach_avoid_polytope, ...
          opt_input_vector_at_vertices, ...
          varargout] = ...
          computeUnderapproxStochReachAvoidSetPwl(...
                                        sys, ...
                                        target_tube, ...
                                        init_safe_set, ...
                                        prob_thresh_of_interest, ...
                                        set_of_direction_vectors,...
                                        varargin)
    method = 'best'; % 'fast'/'best'
    tol_bisection = 1e-2;
    desired_accuracy = 1e-2;
    PSoptions = psoptimset('Display','off');

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

    
    % Get half space representation of the target tube and time horizon
    [concat_target_tube_A, concat_target_tube_b] = target_tube.concat();
    time_horizon = length(target_tube);
    
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

    % Get the no_of_direction vectors
    assert( size(set_of_direction_vectors, 1) == sys.state_dim, ...
            'SReachTools:invalidArgs', ...
            ['set_of_direction_vectors should be a collection of ',...
             ' n-dimensional column vectors.']);
    assert( size(set_of_direction_vectors, 2) >= 2, ...
            'SReachTools:invalidArgs', ...
            'set_of_direction_vectors should at least have two directions.');
    % TODO: Do input handling for uniqueness
    no_of_direction_vectors = size(set_of_direction_vectors, 2);

    % Compute H, mean_X_sans_input_sans_initial_state, cov_X_sans_input, and Z.
    % See @LtiSystem\getConcatMats for the notation.
    % GUARANTEES: Gaussian-perturbed LTI system (sys)
    [H, mean_X_sans_input_sans_initial_state, cov_X_sans_input, Z] = ...
        getHmatMeanCovForXSansInput(sys,zeros(sys.state_dim,1),time_horizon);

    %% Computation of xmax and the associated opt open-loop controller
    % TODO: Translate desired_accuracy to piecewise_count
    piecewise_count = 1e3;
    [invcdf_approx_A, invcdf_approx_b, lower_bound_on_deltai] =...
        computeNormCdfInvOverApprox(piecewise_count);

    %% PWL code preps
    % Compute \sqrt{h_i^\top * \Sigma_X_no_input * h_i}
    sigma_vector = diag(sqrt(diag(concat_target_tube_A *...
        cov_X_sans_input * concat_target_tube_A')));
    % Compute M --- the number of polytopic halfspaces to worry about
    no_linear_constraints = size(concat_target_tube_A,1);

    %% Compute the largest radius possible. This will ensure that W_0(x,U) can
    % be given more importance than the radius in the next optimization problem
    chebyshev_center = init_safe_set.chebyCenter;
    R_upper = chebyshev_center.r;
    % Dual norm of a_i in init_safe_set for chebyshev-centering constraint
    dual_norm_of_init_safe_set_A = sqrt(diag(init_safe_set.A*init_safe_set.A')); 
    %% Solve the optimization problem to compute a point that is deepest inside
    %% the desired set
    % maximize W_0(x,U) - 0.01 * R 
    % subject to
    %   W_0(x,U) >= prob_thresh_of_interest 
    %                      [Computation of W_0(x,U) is done via
    %                       chance-constraint reformulation and risk allocation]
    %   x\in InitState and center of a circle in InitState with radius R
    %   U\in U^N
    %   R>0
    cvx_begin quiet
        variable mean_X(sys.state_dim * time_horizon, 1);
        variable deltai(no_linear_constraints, 1);
        variable norminvover(no_linear_constraints, 1);
        variable xmax(sys.state_dim, 1);
        variable Umax(sys.input_dim * time_horizon,1);
        variable R nonnegative;
    
        maximize ((1-sum(deltai)) + 0.1* (R/R_upper))
        subject to
            % Chance-constraint reformulation of the safety prob (1-6)
            % (1) Slack variables that are bounded below by the piecewise linear
            % approximation of \phi^{-1}(1-\delta)
            for deltai_indx=1:no_linear_constraints
                norminvover(deltai_indx) >= invcdf_approx_A.*...
                    deltai(deltai_indx) + invcdf_approx_b; 
            end
            % (2) Trajectory
            mean_X == Z * xmax + H * Umax +...
                mean_X_sans_input_sans_initial_state;
            % (3) Ono's type of reformulation of chance constraints
            concat_target_tube_A * mean_X + sigma_vector * norminvover...
                <= concat_target_tube_b;
            % (4) Lower bound on delta due to pwl's domain
            deltai >= lower_bound_on_deltai;
            % (5) Upper bound on delta to ensure \phi^{-1}(1-\delta) is concave
            deltai <= 0.5;
            % (6) W_0(x,U) >= alpha
            1-sum(deltai) >= prob_thresh_of_interest
            % Chebyshev-centering constraints
            for i = 1:length(init_safe_set.A)
                init_safe_set.A(i,:) * xmax +...
                    R*dual_norm_of_init_safe_set_A(i)...
                    <= init_safe_set.b(i)
            end
            init_safe_set.Ae * xmax == init_safe_set.be
            % Input constraints
            concat_input_space_A * Umax <= concat_input_space_b;
    cvx_end
    if strcmpi(cvx_status, 'Solved')
        max_underapprox_reach_avoid_prob = 1-sum(deltai);
        opt_input_vector_for_xmax = Umax;
        % Now that we have computed xmax, make sure that the user-provided
        % set_of_direction_vectors will span the affine hull.
        % Compute xmax + directions \forall directions
        check_set_of_direction_vectors =...
            (repmat(xmax,1,no_of_direction_vectors) + set_of_direction_vectors);
        % Check if all of these vectors satisfy the equality constraints of
        % init_safe_set
        are_all_xmaxPlusdirs_within_init_safe_set_eq =...
            all(all(abs(init_safe_set.Ae * check_set_of_direction_vectors -...
                init_safe_set.be)<eps));
        assert( are_all_xmaxPlusdirs_within_init_safe_set_eq,...
            'SReachTools:invalidArgs', ...
            ['set_of_direction_vectors+xmax should lie in the affine hull ',...
             'intersecting safe_set.']);
    else
        % Failed to find a single point within the polytope (if there is at
        % least one x for which W(x,U)= prob_thresh_of_interest, then
        % R would have been >=0, and the above problem would have been feasible)
        max_underapprox_reach_avoid_prob = 0;
        opt_input_vector_for_xmax = nan(sys.input_dim * time_horizon,1);    
    end
    % If non-trivial underapproximative stochastic reach-avoid polytope
    if max_underapprox_reach_avoid_prob < prob_thresh_of_interest
        % Stochastic reach-avoid underapproximation is empty and no admissible
        % open-loop policy exists
        fprintf(['Polytopic underapproximation does not exist for alpha =', ...
                 ' %1.2f since W(x_max) = %1.3f.\n\n'], ...
                 prob_thresh_of_interest, ...
                 max_underapprox_reach_avoid_prob);

        % Assigning the outputs to trivial results
        underapprox_stoch_reach_avoid_polytope = Polyhedron();
        opt_input_vector_at_vertices = nan(...
                                          sys.input_dim * time_horizon, ...
                                          no_of_direction_vectors);        
        xmax = nan(sys.state_dim, 1);
        opt_theta_i = nan(1, no_of_direction_vectors);
        opt_reachAvoid_i = nan(1, no_of_direction_vectors);
        opt_input_vector_at_vertices =...
            nan(sys.input_dim * time_horizon, no_of_direction_vectors);
        vertex_underapprox_polytope = nan(sys.state_dim,...
            no_of_direction_vectors);
        R = nan;
    else
        % Stochastic reach-avoid underapproximation is non-trivial
        fprintf(['Polytopic underapproximation exists for alpha = %1.2f ', ...
                 'since W(x_max) = %1.3f.\n\n'], ...
                 prob_thresh_of_interest, ...
                 max_underapprox_reach_avoid_prob);

        % Pre-allocation of relevant outputs
        opt_input_vector_at_vertices =...
            zeros(sys.input_dim * time_horizon, no_of_direction_vectors);
        opt_theta_i = zeros(1, no_of_direction_vectors);
        opt_reachAvoid_i = zeros(1, no_of_direction_vectors);

        %% Iterate over all direction vectors + xmax
        for direction_index = 1: no_of_direction_vectors
            % Get direction_index-th direction in the hyperplane
            direction = set_of_direction_vectors(:,direction_index);

            %fprintf(['Analyzing direction (shown transposed) ', ...
                     %':%d/%d\n'],direction_index,no_of_direction_vectors);
            %disp(direction');
            fprintf('Analyzing direction :%2d/%2d\n', direction_index, ...
                no_of_direction_vectors);

            %% Solve the optimization problem to compute the boundary point of
            %% the desired underapproximative stochastic reach-avoid set
            % maximize theta
            % subject to
            %   boundary_point = xmax + theta * direction
            %   U\in U^N
            %   boundary_point\in InitState
            %   W_0( boundary_point, U) >= prob_thresh_of_interest 
            %              [As before, the computation of W_0(x,U) is done via
            %               chance-constraint reformulation and risk allocation]
            cvx_begin quiet
                variable mean_X(sys.state_dim * time_horizon, 1);
                variable deltai(no_linear_constraints, 1);
                variable norminvover(no_linear_constraints, 1);
                variable theta nonnegative;
                variable boundary_point(sys.state_dim, 1);
                variable U_vector(sys.input_dim * time_horizon,1);

                maximize theta
                subject to
                    % Define boundary point
                    boundary_point == xmax + theta *direction
                    % Safe boundary point
                    init_safe_set.A * boundary_point <= init_safe_set.b
                    init_safe_set.Ae * boundary_point == init_safe_set.be
                    % Input constraints
                    concat_input_space_A * U_vector <= concat_input_space_b;
                    % Chance-constraint reformulation of safety prob(1-6)
                    % (1) Slack variables that are bounded below by the
                    % piecewise linear approximation of \phi^{-1}(1-\delta)
                    for deltai_indx=1:no_linear_constraints
                        norminvover(deltai_indx) >= invcdf_approx_A.*...
                            deltai(deltai_indx) + invcdf_approx_b; 
                    end
                    % (2) Trajectory
                    mean_X == Z * boundary_point + H * U_vector +...
                        mean_X_sans_input_sans_initial_state;
                    % (3) Ono's type of reformulation of chance constraints
                    concat_target_tube_A * mean_X + sigma_vector* norminvover...
                        <= concat_target_tube_b;
                    % (4) Lower bound on delta due to pwl's domain
                    deltai >= lower_bound_on_deltai;
                    % (5) Upper bound on delta to ensure \phi^{-1}(1-\delta) is
                    % concave
                    deltai <= 0.5;
                    % (6) W_0(x,U) >= alpha
                    1-sum(deltai) >= prob_thresh_of_interest
            cvx_end
            opt_theta_i(direction_index) = theta;
            opt_input_vector_at_vertices(:,direction_index) = U_vector;
            opt_reachAvoid_i(direction_index) = 1 - sum(deltai);
            vertex_underapprox_polytope(:, direction_index)= boundary_point;
        end
        underapprox_stoch_reach_avoid_polytope = Polyhedron('V', ...
            vertex_underapprox_polytope');
    end
    varargout{1} = xmax;
    varargout{2} = opt_input_vector_for_xmax;
    varargout{3} = max_underapprox_reach_avoid_prob;
    varargout{4} = opt_theta_i;
    varargout{5} = opt_reachAvoid_i;
    varargout{6} = R;
    varargout{7} = vertex_underapprox_polytope;
end



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
    %fprintf(['OptRAProb | OptTheta | LB_theta | UB_theta | ',...
             %'OptInp^2 | Exit reason\n']); %10 characters b/n | |

    %% Store opt_input_vector_at_vertices(:,direction_index) for
    %% updates via this bisection
    %opt_inputs_so_far =...
        %opt_input_vector_at_vertices(:,direction_index);
    %mean_X_sans_input = ...
        %Z * boundary_point + mean_X_sans_input_sans_initial_state;
    %opt_reachAvoid_prob_so_far =...
        %computeReachAvoidProb(opt_inputs_so_far, ...
            %mean_X_sans_input,...
            %cov_X_sans_input, ...
            %H, ...
            %concat_target_tube_A, ...
            %concat_target_tube_b, ...
            %desired_accuracy);
    %opt_theta_so_far = theta_lower;
    %% While the gap remains above tolerance_bisection
    %while (theta_upper - theta_lower) > tol_bisection
        %% Set phi as the middle point of the interval
        %phi = (theta_upper+theta_lower)/2;
        %% Candidate initial state provides mean_X_sans_input
        %% See @LtiSystem/getConcatMats for notation
        %candidate_initial_state = xmax + phi * direction;
        %mean_X_sans_input = Z * candidate_initial_state + ...
            %mean_X_sans_input_sans_initial_state;
        %% Compute the maximum reach-avoid prob and the
        %% associated open-loop opt controller starting from
        %% candidate_initial_state using Genz+patternsearch
        %% Add desired_accuracy to prob_thresh_of_interest so that we
        %% are guaranteed to meet the desired threshold 
        %negLogReachAvoidProbGivenInputVectorFeas= @(input_vector)...
            %-log(min(computeReachAvoidProb(input_vector, ...
                        %mean_X_sans_input, ...
                        %cov_X_sans_input, ...
                        %H, ...
                        %concat_target_tube_A, ...
                        %concat_target_tube_b, ...
                        %desired_accuracy),...
                     %prob_thresh_of_interest+desired_accuracy));
        
        %% Compute the opt admissible input_vector that ensures
        %% W_0(x,U)\geq \alpha, given x_0
        %%
        %% minimize
        %% -log(min(ReachAvoidProb(U),prob_thresh_of_interest))
        %% subject to
        %%        concat_input_space_A * U \leq concat_input_space_b
        %[opt_input_at_this_step, patternsearch_opt_value]= ...
              %patternsearch(...
                %negLogReachAvoidProbGivenInputVectorFeas,...
                %opt_inputs_so_far, ...
                %concat_input_space_A, ...
                %concat_input_space_b, ...
                %[],[],[],[],[], ...
                %PSoptions);
        
        %% Recompute the optimal value just to be sure
        %revaluate_opt_prob =...
            %computeReachAvoidProb(opt_input_at_this_step, ...
                %mean_X_sans_input,...
                %cov_X_sans_input, ...
                %H, ...
                %concat_target_tube_A, ...
                %concat_target_tube_b, ...
                %desired_accuracy);

        %if revaluate_opt_prob - prob_thresh_of_interest >=...
                %-tol_bisection && exp(-patternsearch_opt_value) -...
                %prob_thresh_of_interest >= -tol_bisection
            %% Update values only if it is above the thresh
            %opt_reachAvoid_prob_so_far = revaluate_opt_prob;
            %opt_inputs_so_far = opt_input_at_this_step;
            %opt_theta_so_far = phi;
            %exitmessage = ' Feasible';
        %else
            %% Update nothing since Genz+patternsearch failed to
            %% converge
            %exitmessage = sprintf(' Infeasible (%1.3f,%1.3f)',...
                %exp(-patternsearch_opt_value),...
                %revaluate_opt_prob);
        %end
        %% Print current solution
        %fprintf(strcat(...
            %'  %1.4f  |  %1.4f  |  %1.4f  |  %1.4f  |  %1.4f  |  ', ...
            %exitmessage, ' \n'), ...
                   %opt_reachAvoid_prob_so_far, ...
                   %opt_theta_so_far, ...
                   %theta_lower, ...
                   %theta_upper, ...
                   %opt_inputs_so_far'*opt_inputs_so_far);

        %% Update bounds
        %if  strcmpi(exitmessage, ' Feasible')
            %% Theta proved to be feasible
            %theta_lower = phi;
        %else
            %% Shrink the upper bound
            %theta_upper = phi;
        %end
    %end        
    %opt_theta_i(direction_index) = opt_theta_so_far;
    %opt_input_vector_at_vertices(:,direction_index) = opt_inputs_so_far;
    %opt_reachAvoid_i(direction_index) = opt_reachAvoid_prob_so_far;
    %vertex_underapprox_polytope(:, direction_index)= xmax +...
        %opt_theta_so_far*direction;
%end
