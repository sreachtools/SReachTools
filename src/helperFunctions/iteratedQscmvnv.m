function prob = iteratedQscmvnv(qscmvnv_cov, ...
                                qscmvnv_lb, ...
                                qscmvnv_coeff_matrix, ...
                                qscmvnv_ub, ...
                                desired_accuracy, ...
                                warning_iteration)
% SReachTools/helperFunctions/iteratedQscmvnv: Wrapper for Genz's algorithm to
% compute the integral of a Gaussian over an intersection of halfspaces up to a
% desired_accuracy
% ============================================================================
% 
% This function computes the integral of a zero-mean Gaussian \eta with
% covariance matrix qscmvnv_cov over the set of inequalities 
%       qscmvnv_lb <= qscmvnv_coeff_matrix * \eta <= qscmvnv_ub
% Genz's algorithm (qscmvnv) is called iteratively to attain a desired accuracy.
%
% Usage: 
% -----
%
% % To integrate a Gaussian of mean Gauss_mean and covariance Gauss_cov,
% % over a target_set = Polyhedron('H', [A b]), do the following (for an accuracy
% % of 1e-4 and warnings after 10 iterations)
%
% prob=iteratedQscmvnv(Gauss_cov, ...
%                      repmat(-Inf, [size(target_set.A, 1), 1]);
%                      target_set.A;
%                      target_set.b - target_set.A * Gauss_mean;
%                      1e-4, ...
%                      10)
%
% ============================================================================
% 
% prob = iteratedQscmvnv(qscmvnv_cov, ...
%                        qscmvnv_lb, ...
%                        qscmvnv_coeff_matrix, ...
%                        qscmvnv_ub, ...
%                        desired_accuracy, ...
%                        warning_iteration)
%
% Inputs:
% -------
% qscmvnv_cov         - Covariance matrix 
% qscmvnv_lb, qscmvnv_coeff_matrix, qscmvnv_ub 
%                     - The set whose set-membership probability is of interest
%                       qscmvnv_lb <= qscmvnv_coeff_matrix * x 
%                       qscmvnv_ub => qscmvnv_coeff_matrix * x 
% desired_accuracy    - Accuracy of the integral evaluation [If unsure, use 1e-8
%                       if sys.state_dim <= 4, and 1e-3 otherwise]
% warning_iteration   - No. of iterations after which warning should be provided
%                       to make the user aware that the accuracy setting might
%                       be too high [If unsure, use 10]
%
% Outputs:
% --------
%   prob              - Integral of the Gaussian over the set specified
%
% Notes:
% ------
% * NOT ACTIVELY TESTED: TODO
% * NO INPUT HANDLING: For computational speed.
% * MATLAB DEPENDENCY: Uses MATLAB's Statistics and Machine Learning Toolbox.
%                      Need normpdf, norminv, normcdf for Genz's algorithm
% * Uses Genz's algorithm for integral of multivariate Gaussian, qscmvnv, that
%   can be found at
%   http://www.math.wsu.edu/faculty/genz/software/matlab/qscmvnv.m. 
%      * This function has been included in SReachTools/src/helperFunctions. 
%      * This quasi-Monte-Carlo simulations and Cholesky decompostion-based
%        algorithm is driven to provide a desired accuracy by appropriately
%        increasing the number of Monte Carlo particles used.
% * In the event, the integral is below the desired_accuracy, prob is set to
%   desired_accuracy. This is to allow to take log of the prob, if desired.
% * In case, the target set is a hyper-cuboid and the state_dim < 25,
%   then use mvncdf instead.
%
% ============================================================================
%
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
%
%

    %% QSCMVNV in a loop using the error estimate
    error_quadrature = 10;
    points_base = 10;
    points_power = 1;
    while abs(error_quadrature)>desired_accuracy
        [temp_probability, error_quadrature] = ...
                                       qscmvnv(points_base^points_power, ...
                                               qscmvnv_cov, ...
                                               qscmvnv_lb, ...
                                               qscmvnv_coeff_matrix, ...
                                               qscmvnv_ub);
        % Rounding off the integral to the desired accuracy
        temp_probability = round(temp_probability/desired_accuracy) *...
                            desired_accuracy;
        if points_power > warning_iteration
            warning(['Exceeded %d iterations --- Required accuracy: %1.2e |', ...
                     ' Current error: %1.2e\n'], ...
                    warning_iteration, ...
                    desired_accuracy, ...
                    error_quadrature);
        end        
        points_power=points_power+1;
    end
    if points_power > warning_iteration
        warning('Took %d iterations\n',points_power-1)
    end
    % If temp_probability< desired_accuracy, then set reach_avoid_probability to
    % desired_accuracy
    prob = max(temp_probability, desired_accuracy);
end
