function [lb_stoch_reach_avoid, optimal_input_vector, optimal_input_gain, input_satisfaction_prob, risk_alloc] =...
    computeClosedLoopCcLowerBoundStochReachAvoidPwlRisk( ...
        sys,...
        time_horizon,...
        concat_input_space_A, ... 
        concat_input_space_b, ...
        concat_target_tube_A, ... 
        concat_target_tube_b, ...
        H, ...
        G, ...
        mean_X_sans_input, ...
        cov_concat_disturb, ...
        max_state_violation_prob,...
        max_input_violation_prob,...
        desired_accuracy)
% Solve the stochastic reach-avoid problem (lower bound on the probability and
% an open-loop controller synthesis) using chance-constrained convex
% optimization optimization (Internal function --- assumes arguments are all ok)
% =============================================================================
%
% computeCcLowerBoundStochReachAvoidPwlRisk implements the chance-constrained
% convex underapproximation to the terminal hitting-time stochastic reach-avoid
% problem discussed in
%
% K. Lesser, M. Oishi, and R. Erwin, "Stochastic reachability for control of
% spacecraft relative motion," in IEEE Conference on Decision and Control (CDC),
% 2013.
%
% and reformulated in
%
% A. Vinod, V. Sivaramakrishnan, and M. Oishi, CSS-L, 2018 (submitted) TODO
%
% USAGE: This function is intended for internal use as it does not sanitize the
% inputs. Please use getLowerBoundStochReachAvoid instead.
%
% =============================================================================
%   [lb_stoch_reach_avoid, optimal_input_vector] =...
%       computeCcLowerBoundStochReachAvoidPwlRisk( ...
%           sys, ...
%           time_horizon, ...
%           concat_input_space_A, ... 
%           concat_input_space_b, ...
%           concat_target_tube_A, ... 
%           concat_target_tube_b, ...
%           H, ...
%           mean_X_sans_input, ...
%           cov_X_sans_input, ...
%           desired_accuracy)
% 
% Inputs:
% -------
%   sys                   - LtiSystem object
%   time_horizon          - Time horizon (N) with the control provided from 0 to N-1
%   concat_input_space_A, 
%    concat_input_space_b - (A,b) Halfspace representation for the
%                            polytope U^{time_horizon} set.        
%   concat_target_tube_A, 
%    concat_target_tube_b - (A,b) Halfspace representation for the target tube
%                            from t=1 to time_horizon.  For example, the
%                            terminal reach-avoid problem requires a polytope of
%                            the form safe_set^{time_horizon-1} x target_set.        
%   H                     - Concatenated input matrix (see
%                            @LtiSystem/getConcatMats for the notation used)
%   mean_X_sans_input     - Mean of X without the influence of the input
%   cov_X_sans_input      - Covariance of X without the influence of the input
%   desired_accuracy      - Desired accuracy for the optimal stochastic
%                           reach-avoid probability
%
% Outputs:
% --------
%   lb_stoch_reach_avoid - Lower bound on the terminal-hitting stochastic
%                          reach avoid problem computed using Fourier
%                          transform and convex optimization
%   optimal_input_vector - Optimal open-loop policy
%                          ((sys.input_dim) *
%                          time_horizon)-dimensional vector 
%                          U = [u_0; u_1; ...; u_N] (column vector)
%
% See also getLowerBoundStochReachAvoid,
% computeCcLowerBoundStochReachAvoidIterRisk, and
% computeFtLowerBoundStochReachAvoid.
%
% Notes:
% ------
% * NOT ACTIVELY TESTED: Builds on other tested functions.
% * MATLAB DEPENDENCY: Uses MATLAB's Statistics and Machine Learning Toolbox.
%                      Needs norminv
% * EXTERNAL DEPENDENCY: Uses MPT3 and CVX (optional)
%                      Needs MPT3 for defining a controlled system and the
%                      definition of the safe and the target (polytopic) sets
%                      Needs CVX to setup a convex optimization problem that
%                      initializes the patternsearch-based optimization. If CVX
%                      is unavailable, the user may provide a guess for the
%                      initialization.
% * See @LtiSystem/getConcatMats for more information about the
%     notation used.
% 
% ============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 
%

    %% House-keeping for the event where probability <0.5 
    % (TODO: Remove once MILP is encoded)
    lb_stoch_reach_avoid = -1;
    optimal_input_vector = nan(sys.input_dim * time_horizon,1);
    
    %% Compute M --- the number of polytopic halfspaces to worry about
    n_lin_state = size(concat_target_tube_A,1);
    n_lin_input = size(concat_input_space_A,1);
    
    %% Obtain the piecewise linear overapproximation of norminvcdf in the
    %interval [0,0.44]
    [invcdf_approx_m, invcdf_approx_c, lb_deltai,useful_knots] =...
        computeNormCdfInvOverApprox(max(max_state_violation_prob,max_input_violation_prob),...
            desired_accuracy,...
            n_lin_state);
    n_pwa = length(invcdf_approx_m);
    
    d_vector = sdpvar(sys.input_dim * time_horizon, 1);
    mean_X = sdpvar(sys.state_dim * time_horizon, 1);

    deltai = sdpvar(n_lin_state, 1);
    norm_state_replace_slack = sdpvar(n_lin_state, 1);
    slack_cc_sqrt_state = sdpvar(n_lin_state, n_pwa, 'full');
    slack_cc_sqrt_state_iter = sdpvar(n_lin_state, n_pwa, 'full');
    slack_reverse_state = sdpvar(n_lin_state, n_pwa, 'full');

    gammai = sdpvar(n_lin_input, 1);
    norm_input_replace_slack = sdpvar(n_lin_input, 1);
    slack_cc_sqrt_input = sdpvar(n_lin_input, n_pwa, 'full');
    slack_cc_sqrt_input_iter = sdpvar(n_lin_input, n_pwa, 'full');
    slack_reverse_input = sdpvar(n_lin_input, n_pwa, 'full');

    % M_matrix is a decision variable but it is typically buried inside
    % norms
    M_matrix = sdpvar(sys.input_dim * time_horizon, sys.dist_dim * time_horizon, 'full');
    for time_indx = 1:time_horizon - 1
        M_matrix((time_indx-1)*sys.input_dim + 1:time_indx*sys.input_dim, (time_indx-1)*sys.dist_dim+1:end) =0; 
    end

    constraints = [(mean_X == mean_X_sans_input + H * d_vector):'Mean trajectory constraint',...
                   (lb_deltai <= deltai <= max_state_violation_prob):['Risk allocation in [' num2str(lb_deltai) ',' num2str(max_state_violation_prob) '] state'],...
                   (lb_deltai <= gammai <= max_input_violation_prob):['Risk allocation in [' num2str(lb_deltai) ',' num2str(max_input_violation_prob) '] input'],...
                   (sum(deltai) <= max_state_violation_prob):'Max tolerance in state',...
                   (sum(gammai) <= max_input_violation_prob):'Max tolerance in input',...
                   (slack_cc_sqrt_input >= 0):'Positive slack variables CC-sqrt input',...
                   (slack_cc_sqrt_state >= 0):'Positive slack variables CC-sqrt state',...
                   (slack_reverse_state >= 0):'Reverse cvx slack variable in input',... 
                   (slack_reverse_input >= 0):'Reverse cvx slack variable in state'];
    % Reverse convex inequality of slack_cc_sqrt_X^2 - slack_cc_X >= 0 is
    % enforced via slack_cc_X - slack_cc_sqrt_X^2 <= 0 (concave <= 0) by
    % linearizing slack_cc_sqrt_X^2 to get an overapproximation and then
    % making slack_cc_X - overapproximation <= 0. This does not guarantee
    % that slack_cc_X <= slack_cc_sqrt_X^2.


    % Replace all norms with a slack variable whose error is minimized in
    % the optimization problem
    for input_indx = 1:n_lin_input
        constraints = [constraints,...
                       (norm(concat_input_space_A(input_indx,:)* M_matrix * chol(cov_concat_disturb)) <= norm_input_replace_slack(input_indx)):'Slack variable for convexity of norms input'];
    end
    for state_indx = 1:n_lin_state
        constraints = [constraints,...
                       (norm(concat_target_tube_A(state_indx,:)* (H * M_matrix + G) * chol(cov_concat_disturb)) <= norm_state_replace_slack(state_indx)):'Slack variable for convexity of norms state'];
    end                        

    % Enforce the chance constraints as the following two constraints:
    % a^T\mu - b + norm_replace * c <= slack_cc_X                   (a)
    %             slack_cc_sqrt_X^2 <= |m| * norm_replace * delta   (b)
    % Note that the second constraint is a hyperbolic cone constraint
    % which can be reformulated as a second order cone constraint.
    % Specifically, (b) is true iff 
    % || [2*slack_cc_sqrt_X;      ||
    % || norm_replace - |m|*delta]||_2 <= norm_replace + |m|*delta.
    % We enforce slack_cc_X <= slack_cc_sqrt_X^2 (a reverse convex
    % constraint) using DC.
    for approx_indx = 1:length(invcdf_approx_m)
        positive_m_value = abs(invcdf_approx_m(approx_indx));
        positive_c_value = invcdf_approx_c(approx_indx);
        for state_indx = 1:n_lin_state
            constraints = [constraints,...
                           (concat_target_tube_A(state_indx,:) * mean_X   - ...
                                concat_target_tube_b(state_indx) +...
                                norm_state_replace_slack(state_indx) * positive_c_value...
                                - slack_cc_sqrt_state_iter(state_indx,approx_indx).^2 ...
                                - 2 * slack_cc_sqrt_state_iter(state_indx,approx_indx).*(slack_cc_sqrt_state(state_indx,approx_indx) - slack_cc_sqrt_state_iter(state_indx,approx_indx))...
                                    <= slack_reverse_state(state_indx,approx_indx)):['LHS-CC state for each PWA ' num2str(approx_indx)],...
                           (norm([2*slack_cc_sqrt_state(state_indx,approx_indx);
                                  norm_state_replace_slack(state_indx) - positive_m_value * deltai(state_indx)])...
                                        <= norm_state_replace_slack(state_indx) ...
                                            + positive_m_value * deltai(state_indx)):['RHS-CC hyperbolic state for each PWA ' num2str(approx_indx) ', ' num2str(state_indx)]];
        end
        for input_indx = 1:n_lin_input
            constraints = [constraints,...
                           (concat_input_space_A(input_indx,:) * d_vector - ...
                                concat_input_space_b(input_indx) +...
                                norm_input_replace_slack(input_indx) * positive_c_value...
                                - slack_cc_sqrt_input_iter(input_indx,approx_indx).^2 ...
                                - 2 * slack_cc_sqrt_input_iter(input_indx,approx_indx).*(slack_cc_sqrt_input(input_indx,approx_indx) - slack_cc_sqrt_input_iter(input_indx,approx_indx))...
                                    <= slack_reverse_input(input_indx,approx_indx)):['LHS-CC input for each PWA ' num2str(approx_indx)],...
                           (norm([2*slack_cc_sqrt_input(input_indx,approx_indx);
                                  norm_input_replace_slack(input_indx) - positive_m_value * gammai(input_indx)]) <=...
                                        norm_input_replace_slack(input_indx)...
                                            + positive_m_value * gammai(input_indx)):['RHS-CC hyperbolic input for each PWA ' num2str(approx_indx) ', ' num2str(input_indx)]];
        end        
    end

    % Show the constraints
    %constraints;

    % Set some options for YALMIP and solver        
    options = sdpsettings('verbose',0,'solver','gurobi','gurobi.BarIterLimit','1e3','savesolveroutput','1');

    % Minimize the error in the replacement of the norm constraints with
    % slack variables
    objective = sum(norm_input_replace_slack) + sum(norm_state_replace_slack);

    tau = sdpvar(1);
    dc_slack_with_tau = tau * (sum(sum(slack_reverse_state))+sum(sum(slack_reverse_input)));

    t = sdpvar(1);
    ccc_closed_solver = optimizer([constraints,objective + dc_slack_with_tau<=t],...
                                  t,options,...
                                  {slack_cc_sqrt_state_iter, slack_cc_sqrt_input_iter, tau},...
                                  {slack_cc_sqrt_state, slack_cc_sqrt_input, objective, dc_slack_with_tau,...
                                   deltai, gammai, d_vector, M_matrix,...
                                   slack_reverse_state, slack_reverse_input,...
                                   norm_state_replace_slack, norm_input_replace_slack});
    disp('YALMIP has created the optimizer');

    %% Difference of convex penalty approach
    tau_initial = 1;     % Weight for the DC constraint violation
    scaling_tau = 2;     % Multiplying factor for the weight
    tau_max = 1e5;       % Saturation for the weight to avoid numerical issues
    iter_max = 35;       % Max number of DC iterations
    iter_count = 1;      % Counter for the iterations
    dc_conv_tol = 1e-6;  % DC convergence threshold

    % Initializations for DC iterative algorithm
    obj_prev = Inf;      
    dc_slack_with_tau_prev = Inf;
    
    % First solve 
    slack_cc_sqrt_state_iter_curr = 3*ones(n_lin_state, n_pwa);
    slack_cc_sqrt_input_iter_curr = 3*ones(n_lin_input, n_pwa);            
    tau_iter = tau_initial;
    [sol, solver_status] = ccc_closed_solver(slack_cc_sqrt_state_iter_curr, slack_cc_sqrt_input_iter_curr, tau_iter);

    % Iteration status analysis
    obj_curr = sol{3};
    dc_slack_with_tau_curr = sol{4};    
    fprintf([' 0. Status: %d | DC: %1.2e ; Obj: %1.2e | ',...
            'DC slack-total sum state: %1.2e , input: %1.2e\n'],...
            solver_status, dc_slack_with_tau_curr, obj_curr, sum(sum(sol{9})),sum(sum(sol{10})));    
    
    if isnan(obj_curr) || solver_status > 0
        disp('Oops! Initial guess didn''t work');
    else
        %The continue criteria is \leq iter_max AND 
        % NOT OF DC stopping criteria in Lipp and Boyd is met)
        while iter_count <= iter_max && ~(abs((obj_prev + dc_slack_with_tau_prev) - (obj_curr + dc_slack_with_tau_curr)) <= dc_conv_tol)

            % Store previous iterations
            obj_prev = obj_curr;
            dc_slack_with_tau_prev = dc_slack_with_tau_curr;

            % Iteration initialization
            slack_cc_sqrt_state_iter_curr = sol{1};
            slack_cc_sqrt_input_iter_curr = sol{2};            
            tau_iter = min(tau_iter * scaling_tau, tau_max);
            [sol, solver_status] = ccc_closed_solver(slack_cc_sqrt_state_iter_curr, slack_cc_sqrt_input_iter_curr, tau_iter);            

            if isnan(obj_curr) || solver_status > 0
                disp('Converged to a infeasible point! Pausing for diagnostics.');
%                 diag_options = sdpsettings('verbose',0,'solver','gurobi','gurobi.BarIterLimit','1e3','savesolveroutput','1');
%                 ccc_closed_solver_diag = optimizer([constraints,objective + dc_slack_with_tau<=t],...
%                                       t,diag_options,...
%                                       {slack_cc_sqrt_state_iter, slack_cc_sqrt_input_iter, tau},...
%                                       {slack_cc_sqrt_state, slack_cc_sqrt_input, objective, dc_slack_with_tau,...
%                                        deltai, gammai, d_vector, M_matrix,...
%                                        slack_reverse_state, slack_reverse_input,...
%                                        norm_state_replace_slack, norm_input_replace_slack});
%                 [a,b,c,d,e,f]=ccc_closed_solver_diag(slack_cc_sqrt_state_iter_curr,slack_cc_sqrt_input_iter_curr,tau_iter);            
                % where f has the gurobi stuff
                keyboard
            else
                % Iteration status analysis
                obj_curr = sol{3};
                dc_slack_with_tau_curr = sol{4};
                exit_condition = abs(obj_curr + dc_slack_with_tau_curr - obj_prev - dc_slack_with_tau_prev);
                fprintf(['%2d. Status: %d | Total error: %1.2e | tau_iter: %6d | ',...
                         'Curr-Prev --- DC: %1.2e; Obj: %1.2e\n',...
                         'DC slack-total sum state: %1.2e, input: %1.2e\n'],...
                         iter_count, solver_status, exit_condition, tau_iter,... 
                         dc_slack_with_tau_prev - dc_slack_with_tau_curr, obj_prev - obj_curr,...
                         sum(sum(sol{9})),sum(sum(sol{10})));    
                iter_count = iter_count + 1;
            end
        end
    end    
    % Extract and display value
    risk_alloc = sol{5};
    lb_stoch_reach_avoid = 1-sum(sol{5});
    input_satisfaction_prob = 1-sum(sol{6});
    optimal_input_vector = sol{7};
    optimal_input_gain = sol{8};
end
