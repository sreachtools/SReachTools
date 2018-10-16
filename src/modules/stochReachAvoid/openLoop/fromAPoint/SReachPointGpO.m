function [lb_stoch_reach, opt_input_vec] = SReachPointGpO(sys, initial_state,...
    safety_tube, desired_accuracy, PSoptions)
% Solve the stochastic reach-avoid problem (lower bound on the probability and 
% an open-loop controller synthesis) using Fourier transform and convex
% optimization
% =============================================================================
%
% computeFtLowerBoundStochReachAvoid implements the Fourier
% transform-based underapproximation to the terminal hitting-time stochastic
% reach-avoid problem discussed in
%
% A. Vinod and M. Oishi, "Scalable Underapproximation for Stochastic
% Reach-Avoid Problem for High-Dimensional LTI Systems using Fourier
% Transforms," in IEEE Control Systems Letters (L-CSS), 2017.
%
% USAGE: This function is intended for internal use as it does not sanitize the
% inputs. Please use getLowerBoundStochReachAvoid instead.  For the use of
% user-provided initial guess, see getUnderapproxStochReachAvoidSet.
%
% =============================================================================
% [lower_bound_stoch_reach_avoid, optimal_input_vec] = ...
%              computeFtLowerBoundStochReachAvoid(sys, ...
%                                                 time_horizon, ...
%                                                 concat_input_space_A, ... 
%                                                 concat_input_space_b, ...
%                                                 concat_safety_tube_A, ... 
%                                                 concat_safety_tube_b, ...
%                                                 H, ...
%                                                 mean_X_sans_input, ...
%                                                 cov_X_sans_input, ...
%                                                 guess_optimal_input_vec, ...
%                                                 desired_accuracy, ...
%                                                 PSoptions)
% 
% Inputs:
% -------
%   sys                         - LtiSystem object describing the system to be
%                                 verified
%   time_horizon                - Time horizon of the stochastic reach-avoid
%                                 problem
%   concat_input_space_A,       
%    concat_input_space_b       - (A,b) Halfspace representation for the
%                                  polytope U^{time_horizon} set.        
%   concat_safety_tube_A,       
%    concat_safety_tube_b       - (A,b) Halfspace representation for the
%                                 target tube. For example, the terminal
%                                 reach-avoid problem requires a polytope of the
%                                 form safe_set^{time_horizon-1} x target_set.        
%   H                           - Concatenated input matrix (see
%                                 @LtiSystem/getConcatMats for the
%                                 notation used)
%   mean_X_sans_input           - Mean of X
%   cov_X_sans_input            - Covariance of X
%   guess_optimal_input_vec  - User provided initial guess for optimal input
%                                 vector [Use '[]' if unavailable]
%   desired_accuracy            - Accuracy expected for the integral of the
%                                 Gaussian random vector X over the
%                                 concatenated_safety_tube [Use 5e-3 if unsure]
%   PSoptions                   - Options for patternsearch [Use '[]' if unsure]
%
% Outputs:
% --------
%   lower_bound_stoch_reach_avoid - Lower bound on the terminal-hitting 
%                                   stochastic reach avoid problem computed 
%                                   using Fourier transform and convex 
%                                   optimization
%   optimal_input_vec          - Optimal open-loop policy
%                                   ((sys.input_dim) *
%                                   time_horizon)-dimensional vector 
%                                   U = [u_0; u_1; ...; u_N] (column vector)
%
% See also verificationOfCwhDynamics, getLowerBoundStochReachAvoid,
% getUnderapproxStochReachAvoidSet, computeReachAvoidProb
%
% Notes:
% ------
% * NOT ACTIVELY TESTED: Builds on other tested functions.
% * MATLAB DEPENDENCY: Uses MATLAB's Global Optimization Toolbox; Statistics and
%                      Machine Learning Toolbox.
%                      Needs patternsearch for gradient-free optimization
%                      Needs normpdf, normcdf, norminv for Genz's algorithm
% * EXTERNAL DEPENDENCY: Uses CVX (optional)
%                      Needs CVX to setup a convex optimization problem that
%                      initializes the patternsearch-based optimization. If CVX
%                      is unavailable, the user may provide a guess for the
%                      initialization.
% * Uses Genz's algorithm (see in src/helperFunctions) instead of MATLAB's
%   Statistics and Machine Learning Toolbox's mvncdf to compute the integral of
%   the Gaussian over a polytope
% * See @LtiSystem/getConcatMats for more information about the
%   notation used.
% 
% ============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 
%

    % Input handling as well as necessary data creation
    [guess_lb_stoch_reach, guess_opt_input_vec, ~, extra_info] =...
        SReachPointCcO(sys, initial_state, safety_tube, desired_accuracy);
    
    % Get all the other necessary data
    concat_safety_tube_A = extra_info.concat_safety_tube_A;
    concat_safety_tube_b = extra_info.concat_safety_tube_b;
    concat_input_space_A = extra_info.concat_input_space_A;
    concat_input_space_b = extra_info.concat_input_space_b;
    H = extra_info.H;
    mean_X_sans_input = extra_info.mean_X_sans_input;
    cov_X_sans_input = extra_info.cov_X_sans_input;
        
    if guess_lb_stoch_reach < 0
        % Chance constrained approach failed to find an optimal open-loop
        % controller => Try just getting the mean to be as safe as possible
        % 
        % minimize ||s||_1
        % subject to 
        %                               mean_X = mean_X_sans_input + H U
        %  concatenated_safety_tube.A * mean_X <= concatenated_safety_tube.b + s   
        %                                                (Softened reach const.)
        %       concatenated_input_space.A * U <= concatenated_input_space.b
        %
        % Here, mean_X_sans_input is concatenated_A_matrix * x_0 + G_matrix *
        % concatenated_mean_vector_for_the_disturbance (see 
        % getHmatMeanCovForXSansInput)        
        %
        % OVERWRITE the guess_opt_input_vec solution from SReachPointCcO
        cvx_begin quiet
            variable U(sys.input_dim * time_horizon,1)
            variable mean_X(sys.state_dim * time_horizon,1)
            variable slack_vars(length(concat_safety_tube_b))
            minimize norm(slack_vars)
            subject to
                mean_X == mean_X_sans_input + H * U;
                concat_safety_tube_A * mean_X <= concat_safety_tube_b...
                                                    + slack_vars;
                concat_input_space_A * U <= concat_input_space_b;
        cvx_end
        if strcmpi(cvx_status,'Solved') && norm(slack_variable) < 1e-2
            % We have an initial solution to work with
            guess_lb_stoch_reach = 0;
        end
    end

    if guess_lb_stoch_reach < 0
        %% No initial solution to work with
        lb_stoch_reach = -1;
        opt_input_vec = nan(sys.input_dim * time_horizon, 1);
    else
        %% Got an initial solution to work with
        % Construct the reach cost function: 
        % 
        % -log(ReachProb(U)) = -log(\int_{safety_tube} \psi_{\bX}(Z; x_0, U) dZ
        negLogReachProbGivenInputVec= @(input_vec)-log(computeReachAvoidProb(...
            input_vec, mean_X_sans_input, cov_X_sans_input, H, ...
            concat_safety_tube_A, concat_safety_tube_b, desired_accuracy));

        % Compute the optimal admissible input_vec that minimizes
        % negLogReachProbGivenInputVec(input_vec) given x_0 | We use
        % patternsearch because negLogReachProbGivenInputVec can be noisy
        % due to its reliance on Monte Carlo simulation
        %
        % minimize -log(ReachProb(U))
        % subject to
        %        concat_input_space_A * U <= concat_input_space_b
        [opt_input_vec, opt_neg_log_reach_prob]= ...
            patternsearch(negLogReachProbGivenInputVec, guess_opt_input_vec, ...
                concat_input_space_A, concat_input_space_b, [],[],[],[],[],...
                PSoptions);
        
        % Compute the lower bound and the optimal open_loop_control_policy
        lb_stoch_reach = exp(-opt_neg_log_reach_prob);
    end
end