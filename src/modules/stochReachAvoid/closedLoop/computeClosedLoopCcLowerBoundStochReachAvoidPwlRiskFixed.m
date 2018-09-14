function [lb_stoch_reach_avoid, optimal_input_vector, optimal_input_gain, input_satisfaction_prob, risk_alloc] =...
    computeClosedLoopCcLowerBoundStochReachAvoidPwlRiskFixed( ...
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
        max_input_violation_prob)
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

    %% Compute M --- the number of polytopic halfspaces to worry about
    n_lin_state = size(concat_target_tube_A,1);
    n_lin_input = size(concat_input_space_A,1);
    
    risk_alloc_state = max_state_violation_prob/n_lin_state * ones(n_lin_state,1);
    risk_alloc_input = max_input_violation_prob/n_lin_input * ones(n_lin_input,1);
%     risk_alloc_state_mat = diag(norminv(1 - risk_alloc_state));
%     risk_alloc_input_mat = diag(norminv(1 - risk_alloc_input));
    risk_alloc_input_val = norminv(1 - risk_alloc_input(1));
    risk_alloc_state_val = norminv(1 - risk_alloc_state(1));
    
    d_vector = sdpvar(sys.input_dim * time_horizon, 1);
    mean_X = sdpvar(sys.state_dim * time_horizon, 1);
    M_matrix = sdpvar(sys.input_dim * time_horizon, sys.dist_dim * time_horizon, 'full');

    for time_indx = 1:time_horizon - 1
        M_matrix((time_indx-1)*sys.input_dim + 1:time_indx*sys.input_dim, (time_indx-1)*sys.dist_dim+1:end) =0; 
    end
    
    constraints = (mean_X == mean_X_sans_input + H * d_vector):'Mean trajectory constraint';
    for state_indx = 1:n_lin_state
        constraints = [constraints,...
                       (concat_target_tube_A(state_indx,:) * mean_X + risk_alloc_state_val * norm(concat_target_tube_A(state_indx,:)* (H * M_matrix + G) * chol(cov_concat_disturb))  <= concat_target_tube_b(state_indx)):'State CC'];
    end

    for input_indx = 1:n_lin_input
        constraints = [constraints,...
                       (concat_input_space_A(input_indx,:) * d_vector + risk_alloc_input_val * norm(concat_input_space_A(input_indx,:)* M_matrix * chol(cov_concat_disturb))  <= concat_input_space_b(input_indx)):'Input CC'];
    end

    % Set some options for YALMIP and solver
    options = sdpsettings('verbose',0,'solver','gurobi','gurobi.BarIterLimit','1e3');
            
    objective = 0;

    sol = optimize(constraints, objective, options);
    
    if sol.problem > 0
        lb_stoch_reach_avoid = -1;
        input_satisfaction_prob = 0;
        optimal_input_vector = nan(sys.input_dim * time_horizon,1);
        risk_alloc = risk_alloc_state;
        optimal_input_gain = nan(sys.input_dim * time_horizon,sys.dist_dim * time_horizon);
    else
        % Extract and display value
        risk_alloc = risk_alloc_state;
        lb_stoch_reach_avoid = 1-sum(risk_alloc_state);
        input_satisfaction_prob = 1-sum(risk_alloc_input);
        optimal_input_vector = value(d_vector);
        optimal_input_gain = value(M_matrix);
    end
end
