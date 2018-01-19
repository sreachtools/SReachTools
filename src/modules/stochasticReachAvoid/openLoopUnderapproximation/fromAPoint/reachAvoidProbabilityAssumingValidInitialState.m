function [reach_avoid_probability] = ...
   reachAvoidProbabilityAssumingValidInitialState(input_vector,...
                                                  mean_X_sans_input,...
                                                  covariance_X_sans_input,...
                                                  H_matrix,...
                                                  concatenated_target_tube_A,...
                                                  concatenated_target_tube_b,...
                                                  desired_accuracy)
% SReach/stochasticReachAvoid/reachAvoidProbabilityAssumingValidInitialState
% =============================================================================
%
% reachAvoidProbabilityAssumingValidInitialState computes the objective
% function of the Fourier transform-based underapproximation of the terminal
% hitting-time stochastic reach avoid problem as discussed in
%
% A. Vinod and M. Oishi, "Scalable Underapproximation for Stochastic
% Reach-Avoid Problem for High-Dimensional LTI Systems using Fourier
% Transforms," in IEEE Control Systems Letters (L-CSS), 2017.
%
% Specifically, reachAvoidProbabilityAssumingValidInitialState computes the
% integral of the Gaussian random vector (concatenated state vector) X over the
% reach-avoid (polytopic) tube safe_set^{time_horizon-1} x target_set. 
%
% USAGE: See computeFtLowerBoundStochasticReachAvoid.
%
% =============================================================================
%
% [reach_avoid_probability] = ...
%    reachAvoidProbabilityAssumingValidInitialState(input_vector,...
%                                                   mean_X_sans_input,...
%                                                   covariance_X_sans_input,...
%                                                   H_matrix,...
%                                                   concatenated_target_tube_A,...
%                                                   concatenated_target_tube_b,...
%                                                   desired_accuracy)
% 
% Inputs:
% -------
%   input_vector               - Concatenated input vector under investigation
%   mean_X_sans_input          - Mean of (X - H_matrix * input_vector)
%   covariance_X_sans_input    - Covariance matrix of X (Since addition of a
%                                constant to a Gaussian doesn't affect the
%                                covariance matrix)
%   H_matrix                   - 
%   concatenated_target_tube_A - concatenated target tube polyhedral definition
%   concatenated_target_tube_b - concatenated target tube polyhedral definition
%   desired_accuracy           - Accuracy expected for the integral of the
%                                Gaussian random vector X over the
%                                concatenated_target_tube
%
% Outputs:
% --------
%   reach_avoid_probability - Reach-avoid probability attained using the given
%                             input_vector
%
% Notes:
% * NOT ACTIVELY TESTED: TODO
% * NO INPUT HANDLING: For computational speed. To be used via
%   computeFtLowerBoundStochasticReachAvoid
% * MATLAB DEPENDENCY: Uses MATLAB's Statistics and Machine Learning Toolbox.
%                      Need normpdf, norminv, normcdf for Genz's algorithm
% * Uses Genz's algorithm for integral of multivariate Gaussian, qscmvnv, that
%   can be found at
%   http://www.math.wsu.edu/faculty/genz/software/matlab/qscmvnv.m. 
%      * This function has been included in SReach/src/helperFunctions. 
%      * This quasi-Monte-Carlo simulations and Cholesky decompostion-based
%        algorithm is driven to provide a desired accuracy by appropriately
%        increasing the number of Monte Carlo particles used.
% * The integral is rounded off to desired_accuracy provided.
% * In the event, the integral is below the desired_accuracy,
%   reach_avoid_probability is set to desired_accuracy. This is to allow to take
%   log of the reach_avoid_probability.
% * mvncdf may be used if (concatenated_target_tube_A,
%   concatenated_target_tube_b) corresponds to a (sys.state_dimension x
%   time_horizon)-dimensional cuboid and (sys.state_dimension x time_horizon) <
%   25.
% 
% =============================================================================
% 
% This function is part of the Stochastic Optimal Control Toolbox.
% License for the use of this function is given in
%      https://github.com/abyvinod/SReach/blob/master/LICENSE
%
%

    % Construct the mean and covariance of the Gaussian random vector X
    mean_X = mean_X_sans_input + H_matrix * input_vector;
    covariance_X = covariance_X_sans_input;

    % Construct the concatenated target tube polytope for qscmvnv
    qscmvnv_polytope_lower_bound = repmat(-Inf,...
                                      [size(concatenated_target_tube_A, 1), 1]);
    qscmvnv_polytope_coeff_matrix = concatenated_target_tube_A;
    qscmvnv_polytope_upper_bound = concatenated_target_tube_b -...
                                            concatenated_target_tube_A * mean_X;

    %% QSCMVNV in a loop using the error estimate
    error_quadrature = 10;
    points_base = 10;
    points_power = 1;
    warning_iteration = 10;
    while abs(error_quadrature)>desired_accuracy
        [temp_probability, error_quadrature] =...
                                       qscmvnv(points_base^points_power,...
                                               covariance_X,...
                                               qscmvnv_polytope_lower_bound,...
                                               qscmvnv_polytope_coeff_matrix,...
                                               qscmvnv_polytope_upper_bound);
        % Rounding off the integral to the desired accuracy
        temp_probability = round(temp_probability/desired_accuracy) *...
                            desired_accuracy;
        if points_power > warning_iteration
            warning(['Exceeded %d iterations --- Required accuracy: %1.2e |',...
                     ' Current error: %1.2e\n'],...
                    warning_iteration,...
                    desired_accuracy,...
                    error_quadrature);
        end        
        points_power=points_power+1;
    end
    if points_power > warning_iteration
        warning('Took %d iterations\n',points_power-1)
    end
    % If temp_probability< desired_accuracy, then set reach_avoid_probability to
    % desired_accuracy
    reach_avoid_probability = max(temp_probability, desired_accuracy);
end
