function [lower_bound_stoch_reach_avoid, optimal_input_vector] = ...
             computeFtLowerBoundStochReachAvoid(sys, ...
                                                time_horizon, ...
                                                concat_input_space_A, ... 
                                                concat_input_space_b, ...
                                                concat_target_tube_A, ... 
                                                concat_target_tube_b, ...
                                                H, ...
                                                mean_X_sans_input, ...
                                                cov_X_sans_input, ...
                                                guess_optimal_input_vector, ...
                                                desired_accuracy, ...
                                                PSoptions)
% SReachTools/stochasticReachAvoid/computeFtLowerBoundStochReachAvoid: Solve 
% the stochastic reach-avoid problem (lower bound on the probability and an 
% open-loop controller synthesis) using Fourier transform and convex 
% optimization (Internal function --- assumes arguments are all ok)
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
% [lower_bound_stoch_reach_avoid, optimal_input_vector] = ...
%              computeFtLowerBoundStochReachAvoid(sys, ...
%                                                 time_horizon, ...
%                                                 concat_input_space_A, ... 
%                                                 concat_input_space_b, ...
%                                                 concat_target_tube_A, ... 
%                                                 concat_target_tube_b, ...
%                                                 H, ...
%                                                 mean_X_sans_input, ...
%                                                 cov_X_sans_input, ...
%                                                 guess_optimal_input_vector, ...
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
%   concat_target_tube_A,       
%    concat_target_tube_b       - (A,b) Halfspace representation for the
%                                 target tube. For example, the terminal
%                                 reach-avoid problem requires a polytope of the
%                                 form safe_set^{time_horizon-1} x target_set.        
%   H                           - Concatenated input matrix (see
%                                 @LtiSystem/getConcatMats for the
%                                 notation used)
%   mean_X_sans_input           - Mean of X
%   cov_X_sans_input            - Covariance of X
%   guess_optimal_input_vector  - User provided initial guess for optimal input
%                                 vector [Use '[]' if unavailable]
%   desired_accuracy            - Accuracy expected for the integral of the
%                                 Gaussian random vector X over the
%                                 concatenated_target_tube [Use 5e-3 if unsure]
%   PSoptions                   - Options for patternsearch [Use '[]' if unsure]
%
% Outputs:
% --------
%   lower_bound_stoch_reach_avoid - Lower bound on the terminal-hitting 
%                                   stochastic reach avoid problem computed 
%                                   using Fourier transform and convex 
%                                   optimization
%   optimal_input_vector          - Optimal open-loop policy
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

    % MATLAB DEPENDENCY CHECK --- patternsearch
    assert(exist('patternsearch','file')==2, ...
           'SReachTools:setup_error', ...
           'This function needs MATLAB''s Global Optimization Toolbox.');

    %% Construct an initial guess for the input_vector
    if isempty(guess_optimal_input_vector)
        % EXTERNAL DEPENDENCY CHECK --- CVX
        assert(exist('cvx_begin','file')==2, ...
               'SReachTools:setup_error', ...
               ['This function uses CVX if no initial guess (input ', ...
                'policy) for the solver is provided. Please get CVX from ', ...
                'http://cvxr.com.']);
        % Approach 1: Optimize the mean state to be within the reach-avoid tube
        % To account for the case where the mean state can not be held within
        % the reach-avoid tube, we use slack variables which should ideally be
        % zero but if positive then should be as least positive as possible
        length_input_vector = sys.input_dim * time_horizon;
        length_state_vector = sys.state_dim * time_horizon;
        % minimize |s|
        % subject to 
        %                                  X = mean_X_sans_input + \mathscr{H} U
        %     concatenated_target_tube.A * X \leq concatenated_target_tube.b
        %                                          + s    (Softened reach-avoid)
        %     concatenated_input_space.A * U \leq concatenated_input_space.b
        % Here, mean_X_sans_input is concatenated_A_matrix * x_0 + G_matrix *
        % concatenated_mean_vector_for_the_disturbance
        cvx_begin quiet
            variable guess_concatentated_input_vector(length_input_vector)
            variable resulting_X_for_xmax(length_state_vector)
            variable slack_variable(length(concat_target_tube_b))
            minimize norm(slack_variable)
            subject to
                resulting_X_for_xmax == mean_X_sans_input + H * ...
                                            guess_concatentated_input_vector
                concat_target_tube_A * resulting_X_for_xmax <= ...
                                            concat_target_tube_b +...
                                            slack_variable
                concat_input_space_A * guess_concatentated_input_vector...
                                <= concat_input_space_b
        cvx_end
        %TODO: Look at the end of this file for alternative approaches
    else
        guess_concatentated_input_vector = guess_optimal_input_vector;
    end

    if all(isnan(guess_concatentated_input_vector))
        %% No initial solution to work with
        lower_bound_stoch_reach_avoid = 0;
        optimal_input_vector = guess_concatentated_input_vector;
    else
        %% Got an initial solution to work with

        % Construct the reach-avoid cost function: 
        % -log(ReachAvoidProbability(U)) = -log(\int_{ReachAvoidTube}
        % \psi_{\bX}(Z; x_0, U) dZ
        negativeLogReachAvoidProbabilityGivenInputVector = ...
          @(input_vector)-log(computeReachAvoidProb(input_vector, ...
                                                    mean_X_sans_input, ...
                                                    cov_X_sans_input, ...
                                                    H, ...
                                                    concat_target_tube_A, ...
                                                    concat_target_tube_b, ...
                                                    desired_accuracy));

        % Compute the optimal admissible input_vector that minimizes
        % negativeLogReachAvoidProbabilityGivenInputVector(input_vector)
        % minimize -log(ReachAvoidProbability(U))
        % subject to
        %        concat_input_space_A * U \leq concat_input_space_b
        % Given x_0
        [optimal_input_vector, optimal_negative_log_reach_avoid_prob]= ...
              patternsearch(negativeLogReachAvoidProbabilityGivenInputVector,...
                            guess_concatentated_input_vector, ...
                            concat_input_space_A, ...
                            concat_input_space_b, ...
                            [],[],[],[],[], ...
                            PSoptions);
        
        % Compute the lower bound and the optimal open_loop_control_policy
        lower_bound_stoch_reach_avoid = ...
                             exp(-optimal_negative_log_reach_avoid_prob);
    end
end

% TODO: We could have a scaling factor later discounting
%discount_vector = 0.1.^[0:length(concat_target_tube_b)-1];
    %minimize discount_vector*slack_variable

% TODO: Approach 2: Use a convex chance-constrained formulation
% originally proposed by K. Lesser in CDC 2013 paper but modified by A.
% P. Vinod in a future publication.

% Approach 3: Center the trajectory as close as possible
%dual_norm_of_safe_set_A = ...
%    sqrt(diag(concat_target_tube_A*concat_target_tube_A')); 
%cvx_begin quiet
    %variable guess_concatentated_input_vector(length_input_vector)
    %variable resulting_X_for_xmax(length_state_vector)
    %variable R
    %maximize R-0.01*(norm(guess_concatentated_input_vector))
    %subject to
        %resulting_X_for_xmax == mean_X_sans_input + H * ...
                                    %guess_concatentated_input_vector
        %concat_target_tube_A * resulting_X_for_xmax...
                        %+ R * dual_norm_of_safe_set_A <= ...
                                    %concat_target_tube_b
        %concat_input_space_A * ...
         %guess_concatentated_input_vector <= concat_input_space_b
        %R >= 0
%cvx_end
