function prob = getProbReachSet(sys,...
                                initial_state,...
                                target_set,...
                                target_time,...
                                desired_accuracy)
% SReachTools/forwardStochasticReach/getProbReachSet: Compute the probability
% that the state at target_time will lie in the target_set. The starting point
% may be a vector or a RandomVector object
% ============================================================================
%
% This function uses getFSRPDMeanCovariance to compute the forward stochastic
% reach probability density (FSRPD) at a future target time. Next, it evaluates
% the integral of the resulting Gaussian over the user-specified target_set (a
% MPT Polyhedron) using iteratedQscmvnv.
%
% USAGE: See examples/forwardStochasticReachCWH.mlx
%
% ============================================================================
% 
% prob = getProbReachSet(sys,...
%                        initial_state,...
%                        target_set,...
%                        target_time,...
%                        desired_accuracy)
%
% Inputs:
% -------
%   sys              - An object of LtiSystem class 
%   initial_state    - Initial state can be a deterministic n-dimensional vector
%                      or a RandomVector object
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
% See also iteratedQscmvnv, getFSRPDMeanCovariance.
%
% Notes:
% ------
% *  In case, the target set is a hyper-cuboid and the state_dimension < 25,
%    then use mvncdf instead.
% ============================================================================
%
% This function is part of the Stochastic Optimal Control Toolbox.
% License for the use of this function is given in
%      https://github.com/abyvinod/SReachTools/blob/master/LICENSE
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

    % Call Genz's algorithm in an iterative approach to compute the probability
    prob=iteratedQscmvnv(cov_x,...
                         qscmvnv_polytope_lower_bound,...
                         qscmvnv_polytope_coeff_matrix,...
                         qscmvnv_polytope_upper_bound,...
                         desired_accuracy,...
                         10);
end
