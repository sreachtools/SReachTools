function [lb_stoch_reach_avoid, optimal_input_vector, optimal_input_gain, risk_alloc_state, risk_alloc_input] =...
    computeClosedLoopCcLowerBoundStochReachAvoidPwlRiskCVX( ...
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
    

        
%     cov_concat_disturb = kron(eye(time_horizon), ...
%                             sys.dist.parameters.covariance);


    
    %% Difference of convex penalty approach
    tau_initial = 1;     % Weight for the DC constraint violation
    scaling_tau = 2;     % Multiplying factor for the weight
    tau_max = 1e5;       % Saturation for the weight to avoid numerical issues
    iter_max = 100;       % Max number of DC iterations
    iter_count = 0;      % Counter for the iterations
    dc_conv_tol = 1e-4;  % DC convergence threshold

    % Initializations for DC iterative algorithm
    obj_curr = Inf;      
    dc_slack_with_tau_curr = Inf;
    
    % First solve 
    slack_cc_sqrt_state_iter = zeros(n_lin_state, n_pwa);
    slack_cc_sqrt_input_iter = zeros(n_lin_input, n_pwa);            
    tau_iter = tau_initial;

    exit_condition = 1;    
    while exit_condition
        % Store previous iterations
        obj_prev = obj_curr;
        dc_slack_with_tau_prev = dc_slack_with_tau_curr;
            
        % The iteration values are updated at the end of the problem
%         tstart = tic;
        cvx_begin quiet
            variable M_matrix(sys.input_dim * time_horizon, sys.dist_dim * time_horizon);
            variable d_vector(sys.input_dim * time_horizon, 1);
            variable mean_X(sys.state_dim * time_horizon, 1);
            % State chance constraint
            variable deltai(n_lin_state, 1) nonnegative;
            variable norm_state_replace_slack(n_lin_state, 1) nonnegative;
            variable slack_cc_sqrt_state(n_lin_state, n_pwa) nonnegative;
            variable slack_reverse_state(n_lin_state, n_pwa) nonnegative;
            % Input chance constraint
            variable gammai(n_lin_input, 1) nonnegative;
            variable norm_input_replace_slack(n_lin_input, 1) nonnegative;
            variable slack_cc_sqrt_input(n_lin_input, n_pwa) nonnegative;
            variable slack_reverse_input(n_lin_input, n_pwa) nonnegative;
            % Minimize slack variable for the norm replacements (epigraph
            % construction) and also the DC programming-based slack constraints
            minimize (sum(norm_input_replace_slack)+sum(norm_state_replace_slack)...
                        + tau_iter * (sum(sum(slack_reverse_state))+sum(sum(slack_reverse_input))));
            subject to
                % Causality constraints on M_matrix
                for time_indx = 1:time_horizon - 1
                    M_matrix((time_indx-1)*sys.input_dim + 1:time_indx*sys.input_dim, (time_indx-1)*sys.dist_dim+1:end) == 0; 
                end
                % Mean trajectory constraint
                mean_X == mean_X_sans_input + H * d_vector;
                % Risk allocation bounds
                lb_deltai <= deltai <= max_state_violation_prob;
                lb_deltai <= gammai <= max_input_violation_prob;
                % Ensure the thresholds are satisfied
                sum(deltai) <= max_state_violation_prob;
                sum(gammai) <= max_input_violation_prob;
                % Positive slack
                slack_cc_sqrt_state >= slack_cc_sqrt_state_iter/2;
                slack_cc_sqrt_input >= slack_cc_sqrt_input_iter/2;
                slack_reverse_state >= 0;
                slack_reverse_input >= 0;
                % Relaxing the norms with their slack variables
                norms(concat_input_space_A* M_matrix * chol(cov_concat_disturb),2,2) <= norm_input_replace_slack;
                norms(concat_target_tube_A* (H * M_matrix + G) * chol(cov_concat_disturb),2,2) <= norm_state_replace_slack;
                % Enforce CC by introducing a slack variable slac_cc_sqrt_X
                %    a^T\mu - b + norm_replace * c <= slack_cc_sqrt_X^2          (a)
                %                slack_cc_sqrt_X^2 <= |m| * norm_replace * delta (b)
                % Note that (a) is a reverse convex constraint. We enforce it by
                % linearizing the RHS of (a) to obtain
                % a^T\mu - b + norm_replace * c - slack_cc_sqrt_X_iter^2 ...
                %   - 2*slack_cc_sqrt_X_iter(slack_cc_sqrt_X - slack_cc_sqrt_X_iter)
                %   <= slack_reverse_state
                % Note that (b) is a hyperbolic cone constraint which can be
                % reformulated as a second order cone constraint. (b) is true iff 
                % || [2*slack_cc_sqrt_X;      ||
                % || norm_replace - |m|*delta]||_2 <= norm_replace + |m|*delta.
                % Feasibility of these constraints implies feasibility of the
                % problem with relaxed norm constraints
                for approx_indx = 1:length(invcdf_approx_m)
                    positive_m_value = abs(invcdf_approx_m(approx_indx));
                    positive_c_value = invcdf_approx_c(approx_indx);

                    % LHS of the state CC
                    concat_target_tube_A * mean_X - concat_target_tube_b...
                        + norm_state_replace_slack * positive_c_value...
                        - slack_cc_sqrt_state_iter(:,approx_indx).^2 ...
                        - 2 * slack_cc_sqrt_state_iter(:,approx_indx).*(slack_cc_sqrt_state(:,approx_indx) - slack_cc_sqrt_state_iter(:,approx_indx))...
                            <= slack_reverse_state(:,approx_indx);
                    % LHS of the input CC
                    concat_input_space_A * d_vector - concat_input_space_b...
                        + norm_input_replace_slack * positive_c_value...
                        - slack_cc_sqrt_input_iter(:,approx_indx).^2 ...
                        - 2 * slack_cc_sqrt_input_iter(:,approx_indx).*(slack_cc_sqrt_input(:,approx_indx) - slack_cc_sqrt_input_iter(:,approx_indx))...
                            <= slack_reverse_input(:,approx_indx);
                    % RHS of the state CC
                    norms([2*slack_cc_sqrt_state(:,approx_indx),...
                           norm_state_replace_slack - positive_m_value * deltai],2,2) <=...
                                    norm_state_replace_slack + positive_m_value * deltai;
                    % RHS of the input CC
                    norms([2*slack_cc_sqrt_input(:,approx_indx),...
                           norm_input_replace_slack - positive_m_value * gammai],2,2) <=...
                                    norm_input_replace_slack + positive_m_value * gammai;
                end
%         t1=toc(tstart);
        cvx_end
%         t2=toc(tstart);
%         disp(t2);
%         disp(t2-t1);
        
        solver_status = cvx_status;
        % Iteration status analysis
        if strcmpi(cvx_status, 'Solved') || strcmpi(cvx_status, 'Inaccurate/Solved')
            obj_curr = sum(norm_input_replace_slack)+sum(norm_state_replace_slack);
            dc_slack_with_tau_curr = tau_iter * (sum(sum(slack_reverse_state))+sum(sum(slack_reverse_input)));    
                
            if iter_count == 0
                % Exit condition is still 1 since we want to do another
                % iteration
                fprintf(' 0. Status: %s | DC: %1.2e ; Obj: %1.2e\n',...
                        solver_status, dc_slack_with_tau_curr, obj_curr);    
            else
                %The continue criteria is \leq iter_max AND 
                % NOT OF DC stopping criteria in Lipp and Boyd is met)
                exit_condition = iter_count <= iter_max &&...
                    ~(abs((obj_prev + dc_slack_with_tau_prev)...
                           - (obj_curr + dc_slack_with_tau_curr)) <= dc_conv_tol);
                % Iteration status analysis
                fprintf(['%2d. Status: %s\nCurrent slack val: %1.2e | tau_iter: %6d | ',...
                         'Curr-Prev --- DC: %1.2e; Obj: %1.2e\n',...
                         'DC slack-total sum state: %1.2e, input: %1.2e\n',...
                         'DC convergence error: %1.2e\n'],...
                         iter_count, solver_status, obj_curr, tau_iter,... 
                         dc_slack_with_tau_prev - dc_slack_with_tau_curr, obj_prev - obj_curr,...
                         sum(sum(slack_reverse_state)),sum(sum(slack_reverse_input)),...
                         abs((obj_prev + dc_slack_with_tau_prev)...
                           - (obj_curr + dc_slack_with_tau_curr)));    
            end    
        else
            disp('Converged to a infeasible point! Pausing for diagnostics.');
            keyboard
        end
        
        % Next iteration initialization
        slack_cc_sqrt_state_iter = slack_cc_sqrt_state;
        slack_cc_sqrt_input_iter = slack_cc_sqrt_input;            
        tau_iter = min(tau_iter * scaling_tau, tau_max);
        % Increment counter 
        iter_count = iter_count + 1;
    end

    if sum(sum(slack_reverse_state))>= dc_conv_tol || sum(sum(slack_reverse_input)) >= dc_conv_tol
        warning('DC program converged, but slack variables does not respect the chance constraints!');
    end
    % Extract and display value
    lb_stoch_reach_avoid = 1-sum(deltai);
    optimal_input_vector = d_vector;
    optimal_input_gain = M_matrix;
    risk_alloc_state = deltai;
    risk_alloc_input = gammai;
end
