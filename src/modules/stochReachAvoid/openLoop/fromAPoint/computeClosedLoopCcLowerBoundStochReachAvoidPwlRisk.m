function [lb_stoch_reach_avoid, optimal_input_vector, optimal_input_gain] =...
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
        cov_X_sans_input, ...
        cov_concat_disturb, ...
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
        computeNormCdfInvOverApprox(0.44,...
            5e-3,...
            n_lin_state);

    %% Solve the feasibility problem
%     cvx_begin quiet
%         variable U_vector(sys.input_dim * time_horizon,1);
%         variable mean_X(sys.state_dim * time_horizon, 1);
%         variable deltai(n_lin_const, 1);
%         variable norminvover(n_lin_const, 1);
%         minimize sum(deltai)
%         subject to
%             mean_X == mean_X_sans_input + H * U_vector;
%             concat_input_space_A * U_vector <= concat_input_space_b;
%             for deltai_indx=1:n_lin_const
%                 norminvover(deltai_indx) >= invcdf_approx_m.*...
%                     deltai(deltai_indx) + invcdf_approx_c; 
%             end
%             concat_target_tube_A * mean_X + sigma_vector * norminvover...
%                 <= concat_target_tube_b;
%             deltai >= lb_deltai;
%             deltai <= 0.5;
%             sum(deltai) <= 0.5;
%     cvx_end
    Delta_input = 1;
    Delta_state = 0.2;
    d_vector = sdpvar(sys.input_dim * time_horizon, 1);
    mean_X = sdpvar(sys.state_dim * time_horizon, 1);
    
    n_pwa = length(invcdf_approx_m);
    
    deltai = sdpvar(n_lin_state, 1);
    sigma_vector_state = sdpvar(n_lin_state, 1);
    sigma_vector_state_slack = sdpvar(n_lin_state, 1);
    slack_cc_state = sdpvar(n_lin_state, n_pwa, 'full');
    
    gammai = sdpvar(n_lin_input, 1);
    sigma_vector_input = sdpvar(n_lin_input, 1);
    sigma_vector_input_slack = sdpvar(n_lin_input, 1);
    slack_cc_input = sdpvar(n_lin_input, n_pwa, 'full');
    
    % M_matrix is a decision variable but it is typically buried inside
    % norms
    M_matrix = sdpvar(sys.input_dim * time_horizon, sys.dist_dim * time_horizon, 'full');
    sqrt_cov_state_matrix = sdpvar(sys.state_dim * time_horizon, sys.dist_dim * time_horizon);
    sqrt_cov_state_matrix = (H * M_matrix + G) * chol(cov_concat_disturb);
    sqrt_cov_input_matrix = sdpvar(sys.input_dim * time_horizon, sys.dist_dim * time_horizon);
    sqrt_cov_input_matrix = M_matrix * chol(cov_concat_disturb);
    
    objective = sum(sigma_vector_input_slack) + sum(sigma_vector_state_slack);%sum(deltai);%sum(sigma_vector_state) + sum(sigma_vector_input);

    for time_indx = 1:time_horizon - 1
        M_matrix((time_indx-1)*sys.input_dim + 1:time_indx*sys.input_dim, (time_indx)*sys.dist_dim+1:end) =0; 
    end
    
    constraints = [(mean_X == mean_X_sans_input + H * d_vector):'Mean trajectory constraint',...
                   (0 <= slack_cc_input):'Positive slack variables input',...
                   (0 <= slack_cc_state):'Positive slack variables state',...
                   (sigma_vector_input <= sigma_vector_input_slack):'Slack variable for convexity of norms input'
                   (sigma_vector_state <= sigma_vector_state_slack):'Slack variable for convexity of norms state'
                   (lb_deltai <= gammai <= 0.43):'Risk allocation in [lb_deltai,0.43] input',...  %0.43 since we can't go all the way to 0.5
                   (lb_deltai <= deltai <= 0.43):'Risk allocation in [lb_deltai,0.43] state',...
                   (sum(deltai) <= Delta_state):'Max tolerance in state',...   %TODO: Is it 0.5 or 1 or necessary?
                   (sum(gammai) <= Delta_input):'Max tolerance in input'];
    
   for approx_indx = 1:length(invcdf_approx_m)
        positive_m_value = abs(invcdf_approx_m(approx_indx));
        positive_c_value = invcdf_approx_c(approx_indx);
        constraints = [constraints,...
                       (concat_input_space_A * d_vector - concat_input_space_b + sigma_vector_input * positive_c_value <= slack_cc_input(:,approx_indx)):['LHS-CC input for each PWA ' num2str(approx_indx)],...
                       (concat_target_tube_A * mean_X   - concat_target_tube_b + sigma_vector_state * positive_c_value <= slack_cc_state(:,approx_indx)):['LHS-CC state for each PWA ' num2str(approx_indx)]];
        for state_indx = 1:n_lin_state
            sigma_vector_state(state_indx) = norm(concat_target_tube_A(state_indx,:)*sqrt_cov_state_matrix);
            constraints = [constraints,...
                           (slack_cc_state(state_indx,approx_indx) <= positive_m_value * sigma_vector_state_slack(state_indx) * deltai(state_indx)):['RHS-CC state for each PWA ' num2str(approx_indx) ', ' num2str(state_indx)]];                       
        end
        
        for input_indx = 1:n_lin_input
            sigma_vector_input(input_indx) = norm(concat_input_space_A(input_indx,:)*sqrt_cov_input_matrix);
            constraints = [constraints,...
                           (slack_cc_input(input_indx,approx_indx) <= positive_m_value * sigma_vector_input_slack(input_indx) * gammai(input_indx)):['RHS-CC input for each PWA ' num2str(approx_indx) ', ' num2str(input_indx)]];                       
        end
    end
               
    % Set some options for YALMIP and solver
    options = sdpsettings('verbose',2);
    
    constraints
    
    % Solve the problem
    sol = optimize(constraints,objective,options);

    % Analyze error flags
    if sol.problem == 0
        % Extract and display value
        lb_stoch_reach_avoid = 1-sum(value(deltai));
        optimal_input_vector = value(d_vector);
        optimal_input_gain = value(M_matrix);
    else
        sol
    end
end
