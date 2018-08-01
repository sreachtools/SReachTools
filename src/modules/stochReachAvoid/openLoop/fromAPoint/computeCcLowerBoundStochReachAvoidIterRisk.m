function [lb_stoch_reach_avoid, optimal_input_vector] =...
    computeCcLowerBoundStochReachAvoidIterRisk( ...
        sys,...
        time_horizon,...
        concat_input_space_A, ... 
        concat_input_space_b, ...
        concat_target_tube_A, ... 
        concat_target_tube_b, ...
        H, ...
        mean_X_sans_input, ...
        cov_X_sans_input, ...
        desired_accuracy)
% Solve the stochastic reach-avoid problem (lower bound on the probability and
% an open-loop controller synthesis) using chance-constrained convex
% optimization optimization (Internal function --- assumes arguments are all ok)
% =============================================================================
%
% computeCcLowerBoundStochReachAvoidIterRisk implements the chance-constrained
% convex underapproximation to the terminal hitting-time stochastic reach-avoid
% problem discussed in
%
% K. Lesser, M. Oishi, and R. Erwin, "Stochastic reachability for control of
% spacecraft relative motion," in IEEE Conference on Decision and Control (CDC),
% 2013.
%
% via 
%
% M. Ono and B. Williams, "Iterative Risk Allocation: A new approach to robust
% Model Predictive Control with a joint chance constraint", in IEEE Conference
% on Decision and Control (CDC), 2008.
%
% USAGE: This function is intended for internal use as it does not sanitize the
% inputs. Please use getLowerBoundStochReachAvoid instead.
%
% =============================================================================
%   [lb_stoch_reach_avoid, optimal_input_vector] =...
%       computeCcLowerBoundStochReachAvoidIterRisk( ...
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
%   sys                  - LtiSystem object describing the system to be verified
%   time_horizon         - Time horizon of the stochastic reach-avoid problem
%   concat_input_space_A,
%    concat_input_space_b- (A,b) Halfspace representation for the
%                           polytope U^{time_horizon} set.        
%   concat_target_tube_A,
%    concat_target_tube_b- (A,b) Halfspace representation for the target tube.
%                           For example, the terminal reach-avoid problem
%                           requires a polytope of the form
%                           safe_set^{time_horizon-1} x target_set.        
%   H                    - Concatenated input matrix (see
%                           @LtiSystem/getConcatMats for the notation used)
%   mean_X_sans_input    - Mean of X without the influence of the input
%   cov_X_sans_input     - Covariance of X without the influence of the input
%   desired_accuracy     - Desired accuracy for the optimal stochastic
%                           reach-avoid probability [If unsure, use 1e-3]
%
% Outputs:
% --------
%   lb_stoch_reach_avoid - Lower bound on the terminal-hitting stochastic
%                          reach avoid problem computed using Fourier
%                          transform and convex optimization
%   optimal_input_vector - Optimal open-loop policy
%                          ((sys.input_dim) * time_horizon)-dimensional 
%                          vector U = [u_0; u_1; ...; u_N] (column vector)
%
% See also verificationOfCwhDynamics,
% getFtLowerBoundStochasticReachAvoid,
% getFtBasedUnderapproxStochReachAvoidSet and
% reachAvoidProbAssumingValidInitialState
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

    %% Ono's converge alpha value
    alpha_on_iter = @(n) 0.7 * (0.98)^n;

    %% House-keeping for the event where the bisection fails to take off
    lb_stoch_reach_avoid = -1;
    optimal_input_vector = nan(sys.input_dim * time_horizon,1);
    
    %% Define when a slack variable will be set to zero
    myeps = 1e-8;

    %% Compute M --- the number of polytopic halfspaces to worry about
    no_linear_constraints = size(concat_target_tube_A,1);

    %% Compute \sqrt{h_i^\top * \Sigma_X_no_input * h_i}
    sigma_vector = sqrt(diag(concat_target_tube_A *...
                                     cov_X_sans_input * concat_target_tube_A'));

    %% Bisection over Delta between 0 and 0.5    
    Delta_lb = 0;
    Delta_ub = 0.5;
    while (Delta_ub - Delta_lb) > desired_accuracy
        Delta = (Delta_ub + Delta_lb)/2;

        %% Prepration for iterated risk allocation
        % Variables (re)initialized
        iter_count = 0;                 % No. of iterations done in Ono's risk
                                        % allocation algorithm
        opt_value_prev = 10000;         % Previous optimal value --- exit cond.
        opt_value = 0;                  % Optimal value --- |slack variables|_1
        % Given Delta, construct delta_i as Delta/M
        delta_vec = Delta/no_linear_constraints * ones(no_linear_constraints,1);

        % Converge to one more decimal precision
        while abs(opt_value - opt_value_prev) >= desired_accuracy/10
            %% Store the previous optimal value
            opt_value_prev = opt_value;

            %% Solve the feasibility problem
            % Construct the back-off (Hessem's term) in the constraints
            scaled_norminv=sigma_vector.*...
                              norminv(ones(no_linear_constraints,1)- delta_vec);
            cvx_begin quiet
                variable U_vector(sys.input_dim * time_horizon,1);
                variable mean_X(sys.state_dim * time_horizon, 1);
                variable slack_variables(no_linear_constraints, 1);
                minimize norm(slack_variables,1);
                subject to
                    mean_X == mean_X_sans_input + H * U_vector;
                    concat_input_space_A * U_vector <= concat_input_space_b;
                    concat_target_tube_A * mean_X + scaled_norminv <= ...
                                         concat_target_tube_b + slack_variables; 
                    slack_variables >= 0;
            cvx_end
            opt_value = norm(slack_variables,1);

            %% Number of active/infeasible constraints via complementary
            %% slackness --- non-zero slack variables imply infeasible \delta_i
            N_active=nnz(slack_variables >= myeps);

            %% Break if N_active is zero or all are active
            if N_active == 0 || N_active == no_linear_constraints
                break;
            end

            %% Compute \alpha value --- decides on how quickly inactive \delta_i
            %% are tightened: 
            alpha = alpha_on_iter(iter_count); 

            %% For all the inactive constraints, \delta_i^+ \gets \alpha
            %% \delta_i + (1-\alpha) (1-normcdf(g-hx/sigma_vec(i)))
            inactive_indx = find(slack_variables < myeps);
            % Compute relevant g-hx^\ast
            correction_deltas = concat_target_tube_b(inactive_indx) - ...
                concat_target_tube_A(inactive_indx,:) * mean_X;
            % Update inactive delta
            delta_vec(inactive_indx) = alpha * delta_vec(inactive_indx) +...
                (1 - alpha) * (1 - normcdf(...
                    correction_deltas./sigma_vector(inactive_indx)));

            %% Collect the headroom gained: \delta_residual \gets \Delta -
            %% \sum_i \delta_i
            delta_residual = Delta - sum(delta_vec);

            %% For all the active constraints, distribute the remaining headroom
            active_indx = find(slack_variables >= myeps);
            % \delta_i \gets \delta_i + \delta_residual/N_active 
            delta_vec(active_indx) = delta_vec(active_indx) +...
                delta_residual / N_active;

            %% Update the iteration count
            iter_count = iter_count + 1;
        end
        fprintf('Done with Delta: %1.4f, N_active: %2d ',...
                Delta,...
                N_active);
        if N_active == 0
            % If no constraint is active, then a feasible input policy has been
            % found --- Dream higher (decrease delta)
            Delta_ub = Delta;
            lb_stoch_reach_avoid = 1-Delta;
            optimal_input_vector = U_vector;
            disp(' Aim for more probability');
        else
            % All constraints are active, then no solution could be found ---
            % Dream lower (increase delta)
            Delta_lb = Delta;
            disp(' Aim for less probability');
        end
    end   
end
