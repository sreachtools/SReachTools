function prob = getProbReachSet(sys,...
                                initial_state,...
                                target_set,...
                                target_time,...
                                desired_accuracy)
% SReach/forwardStochasticReach/getProbReachSet: Compute the probability that
% the state at target_time will lie in the target_set
% ============================================================================
%
% This function is a wrapper for Genz's algorithm to compute the integral of a
% Gaussian over a polytope. Specifically, it uses getFSRPDMeanCovariance to
% compute the forward stochastic reach probability density (FSRPD) at a future
% target time and the evaluates the integral of the resulting Gaussian over the
% user-specified target_set (a MPT Polyhedron).
%
% USAGE: 
%
% ============================================================================
% 
% prob = getProbReachSet(sys,...
%                        initial_state,...
%                        target_set,...
%                        target_time,...
%                        desired_accuracy)
% Inputs:
% -------
%   sys              - An object of LtiSystem class 
%   initial_state    - x_0
%   target_time      - Time of interest (positive scalar)
%   target_set       - Target set 
%   desired_accuracy - Accuracy of the integral evaluation 
%                      [If unsure, use 1e-8 if sys.state_dimension <= 4
%                                      1e-3 otherwise]
%
% Outputs:
% --------
%   prob             - Probability that x_{target_time} lies in target_set
%
% Notes:
% ------
% *  In case, the target set is a hyper-cuboid and the state_dimension < 25,
%    then use mvncdf instead.
% ============================================================================
%
% This function is part of the Stochastic Optimal Control Toolbox.
% License for the use of this function is given in
%      https://github.com/abyvinod/SReach/blob/master/LICENSE
%
%
    [mean_x, cov_x] = getFSRPDMeanCovariance(sys,...
                                             initial_state,...
                                             target_time);
    
    % Construct the concatenated target tube polytope for qscmvnv
    qscmvnv_polytope_lower_bound = repmat(-Inf,...
                                      [size(target_set.A, 1), 1]);
    qscmvnv_polytope_coeff_matrix = target_set.A;
    qscmvnv_polytope_upper_bound = target_set.b - target_set.A * mean_x;

    %% QSCMVNV in a loop using the error estimate
    error_quadrature = 10;
    points_base = 10;
    points_power = 1;
    warning_iteration = 10;
    while abs(error_quadrature)>desired_accuracy
        [temp_probability, error_quadrature] =...
                                       qscmvnv(points_base^points_power,...
                                               cov_x,...
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
    prob = max(temp_probability, desired_accuracy);
end
