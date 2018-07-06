function [underapprox_stoch_reach_avoid_polytope, ...
          optimal_input_vector_at_boundary_points, ...
          varargout] = ...
          computeUnderapproxStochReachAvoidSetPwl(...
                                        sys, ...
                                        target_tube, ...
                                        init_safe_set, ...
                                        probability_threshold_of_interest, ...
                                        no_of_direction_vectors, ...
                                        affine_hull_of_interest_2D)

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

    % Check probability_threshold_of_interest is a scalar in [0.5,1]
    assert( isscalar(probability_threshold_of_interest) &&...
                probability_threshold_of_interest >= 0.5 &&...
                probability_threshold_of_interest <= 1, ...
            'SReachTools:invalidArgs', ...
            'probability_threshold_of_interest must be a scalar in [0.5,1]');

    % Get the set of direction vectors
    % GUARANTEES: sanitized no_of_direction_vectors, affine_hull_of_interest_2D
    set_of_direction_vectors = computeDirectionVectors( ...
                                                 no_of_direction_vectors, ...
                                                 sys.state_dim, ...
                                                 affine_hull_of_interest_2D);

    % Compute H, mean_X_sans_input_sans_initial_state, cov_X_sans_input, and Z.
    % See @LtiSystem\getConcatMats for the notation.
    % GUARANTEES: Gaussian-perturbed LTI system (sys)
    [H, mean_X_sans_input_sans_initial_state, cov_X_sans_input, Z] = ...
        getHmatMeanCovForXSansInput(sys,zeros(sys.state_dim,1),time_horizon);

    %% Computation of xmax and the associated optimal open-loop controller
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
    % Chebyshev-centering constraint
    dual_norm_of_init_safe_set_A = sqrt(diag(init_safe_set.A*init_safe_set.A')); 

    %% Initializing variables
    max_underapprox_reach_avoid_prob = 0;
    optimal_input_vector_for_xmax = nan(sys.input_dim * time_horizon,1);    

    %% Solve the maximize W_0(x,U) with x\in InitState
    R_lower = 0;
    cvx_begin quiet
        variable y(sys.state_dim)
        variable R_upper

        maximize R_upper

        subject to
            R_upper >= 0
            for i = 1:length(init_safe_set.A)
                init_safe_set.A(i,:) * y +...
                    R_upper*dual_norm_of_init_safe_set_A(i)...
                    <= init_safe_set.b(i)
            end
    cvx_end
    disp(R_upper)
    R_tol = 1e-3;
    while abs(R_upper - R_lower) >= R_tol
        R = (R_upper + R_lower)/2;
        cvx_begin quiet
            variable mean_X(sys.state_dim * time_horizon, 1);
            variable deltai(no_linear_constraints, 1);
            variable norminvover(no_linear_constraints, 1);
            variable xmax(sys.state_dim, 1);
            variable Umax(sys.input_dim * time_horizon,1);
    
            minimize sum(deltai)
            subject to
                % Trajectory
                mean_X == Z * xmax + H * Umax +...
                    mean_X_sans_input_sans_initial_state;
                % Input constraints
                concat_input_space_A * Umax <= concat_input_space_b;
                % Slack variables that are bounded below by the piecewise linear
                % approximation of \phi^{-1}(1-\delta)
                for deltai_indx=1:no_linear_constraints
                    norminvover(deltai_indx) >= invcdf_approx_A.*...
                        deltai(deltai_indx) + invcdf_approx_b; 
                end
                % Ono's linear-like constraint 
                concat_target_tube_A * mean_X + sigma_vector * norminvover...
                    <= concat_target_tube_b;
                deltai >= lower_bound_on_deltai;
                deltai <= 0.5;
                %init_safe_set.A * xmax <= init_safe_set.b
                for i = 1:length(init_safe_set.A)
                    init_safe_set.A(i,:) * xmax...
                        + R * dual_norm_of_init_safe_set_A(i) <=...
                            init_safe_set.b(i)
                end
        cvx_end
        if strcmpi(cvx_status, 'Solved')
            % Aim for larger R
            R_lower = R;
            max_underapprox_reach_avoid_prob = 1-sum(deltai);
            optimal_input_vector_for_xmax = Umax;    
        else
            % Aim for smaller R
            R_upper = R;
        end
        disp(R)
    end
    % If non-trivial underapproximative stochastic reach-avoid polytope
    if max_underapprox_reach_avoid_prob < probability_threshold_of_interest
        % Stochastic reach-avoid underapproximation is empty and no admissible
        % open-loop policy exists
        fprintf(['Polytopic underapproximation does not exist for alpha =', ...
                 ' %1.2f since W(x_max) = %1.3f.\n\n'], ...
                 probability_threshold_of_interest, ...
                 max_underapprox_reach_avoid_prob);

        % Assigning the outputs to trivial results
        underapprox_stoch_reach_avoid_polytope = Polyhedron();
        optimal_input_vector_at_boundary_points = nan(...
                                          sys.input_dim * time_horizon, ...
                                          no_of_direction_vectors);        
        optimal_theta_i = zeros(1, no_of_direction_vectors);
        optimal_reachAvoid_i = zeros(1, no_of_direction_vectors);
        optimal_inputs_i = nan(sys.input_dim * time_horizon, ...
                                 no_of_direction_vectors);
        xmax = nan(sys.state_dim, 1);
    else
        % Stochastic reach-avoid underapproximation is non-trivial
        fprintf(['Polytopic underapproximation exists for alpha = %1.2f ', ...
                 'since W(x_max) = %1.3f.\n\n'], ...
                 probability_threshold_of_interest, ...
                 max_underapprox_reach_avoid_prob);

        % For storing boundary points
        optimal_theta_i = zeros(1, no_of_direction_vectors);
        optimal_reachAvoid_i = zeros(1, no_of_direction_vectors);
        optimal_inputs_i = zeros(sys.input_dim * time_horizon, ...
                                 no_of_direction_vectors);

        %% Iterate over all direction vectors + xmax
        for direction_index = 1: no_of_direction_vectors
            % Get direction_index-th direction in the hyperplane
            direction = set_of_direction_vectors(:,direction_index);

            fprintf(['Analyzing direction (shown transposed) ', ...
                     ':%d/%d\n'],direction_index,no_of_direction_vectors);
            disp(direction');

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
                    % Trajectory
                    mean_X == Z * boundary_point + H * U_vector +...
                        mean_X_sans_input_sans_initial_state;
                    % Input constraints
                    concat_input_space_A * U_vector <= concat_input_space_b;
                    % Slack variables that are bounded below by the piecewise
                    % linear approximation of \phi^{-1}(1-\delta)
                    for deltai_indx=1:no_linear_constraints
                        norminvover(deltai_indx) >= invcdf_approx_A.*...
                            deltai(deltai_indx) + invcdf_approx_b; 
                    end
                    % Ono's linear-like constraint 
                    concat_target_tube_A * mean_X + sigma_vector * norminvover...
                        <= concat_target_tube_b;
                    deltai >= lower_bound_on_deltai;
                    deltai <= 0.5;
                    1 - sum(deltai) >= probability_threshold_of_interest
            cvx_end
            optimal_theta_i(direction_index) = theta;
            optimal_inputs_i(:,direction_index) = U_vector;
            optimal_reachAvoid_i(direction_index) = 1 - sum(deltai);
        end
        vertex_underapprox_polytope = xmax +...
                                      optimal_theta_i.*set_of_direction_vectors;
        underapprox_stoch_reach_avoid_polytope = Polyhedron('V', ...
                                                  vertex_underapprox_polytope');
        % Assignment to the respective outputs of this function
        optimal_input_vector_at_boundary_points = optimal_inputs_i;
    end
    varargout{1} = xmax;
    varargout{2} = optimal_input_vector_for_xmax;
    varargout{3} = max_underapprox_reach_avoid_prob;
    varargout{4} = optimal_theta_i;
    varargout{5} = optimal_reachAvoid_i;
end
