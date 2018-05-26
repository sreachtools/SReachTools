function [mean_x, cov_x] = getFSRPDMeanCovariance(sys,...
                                                  initial_state,...
                                                  target_time)
% SReachTools/forwardStochasticReach/getFSRPDMeanCovariance: Compute the mean 
% and the covariance of the state at a time instant in future
% ============================================================================
%
% Computes the mean and the covariance of a Gaussian-perturbed LTI uncontrolled
% system. This function implements Proposition 1 of
%
% A. Vinod, B. HomChaudhuri, and M. Oishi, "Forward Stochastic Reachability
% Analysis for Uncontrolled Linear Systems using Fourier Transforms", In
% Proceedings of the 20th International Conference on Hybrid Systems:
% Computation and Control (HSCC), 2017.
%
% USAGE: See getProbReachSet.m, examples/forwardStochasticReachCWH.mlx.
%
% ============================================================================
% 
% [mean_x, cov_x] = getFSRPDMeanCovariance(sys,...
%                                          initial_state,...
%                                          target_time)
% Inputs:
% -------
%   sys           - An object of LtiSystem class 
%   initial_state - Initial state can be a deterministic n-dimensional vector
%                   or a RandomVector object
%   target_time   - Time of interest (positive scalar)
%
% Outputs:
% --------
%   mean_x        - Mean of the stochastic disturbance
%   cov_x         - Covariance of the stochastic disturbance
%
% See also getProbReachSet
%
% ============================================================================
%
% This function is part of the Stochastic Optimal Control Toolbox.
% License for the use of this function is given in
%      https://github.com/abyvinod/SReachTools/blob/master/LICENSE
%
%

    %% Input handling
    if ~isa(initial_state,'RandomVector')
        % Ensure that initial state is a column vector of appropriate dimension
        assert( size(initial_state,1) == sys.state_dimension &&...
                size(initial_state,2) == 1,...
               'SReachTools:invalidArgs',...
               ['Expected a sys.state_dimension-dimensional column-vector ',...
               'for initial state']);
    end

    % Ensure that the given system has a Gaussian disturbance
    assert(strcmp(class(sys.disturbance),'StochasticDisturbance') &&...
           strcmp(sys.disturbance.type,'Gaussian'),...
          'SReachTools:invalidArgs',...
          'getFSRPDMeanCovariance accepts only Gaussian-perturbed LTI systems');

    % Ensure that the given system is uncontrolled
    assert(sys.input_dimension == 0,...
          'SReachTools:invalidArgs',...
          'getFSRPDMeanCovariance accepts only uncontrolled LTI systems');

    % Ensure that the target_time is a positive scalar
    assert( isscalar(target_time) && target_time > 0,...
           'SReachTools:invalidArgs',...
           'Expected a scalar positive target_time');

    % IID assumption allows to compute the mean and covariance of the
    % concatenated disturbance vector W
    mean_concat_disturb = kron(ones(target_time,1), ...
                               sys.disturbance.parameters.mean);
    cov_concat_disturb =kron(eye(target_time), ...
                             sys.disturbance.parameters.covariance);

    % Controllability matrix for disturbance
    flipped_ctrb_mat_disturb = sys.disturbance_matrix;
    for time_index=1:target_time-1
        % Prepend A times (A^(time_index-1) * F) to the existing 
        % flipped_controllability_matrix_disturbance
        flipped_ctrb_mat_disturb = ...
            [sys.state_matrix * flipped_ctrb_mat_disturb(:,...
                                            1:sys.disturbance_dimension),...
             flipped_ctrb_mat_disturb];
    end

    A_power_target_time = sys.state_matrix^target_time;
    if isa(initial_state,'RandomVector')
        % mu_x = A^\tau * mu_{x_0} + C * mu_{W};
        mean_x = A_power_target_time * initial_state.parameters.mean + ...
                                flipped_ctrb_mat_disturb * mean_concat_disturb;
        % cov_x = A^\tau * cov_{x_0} * (A^\tau)' + C * cov_{W} * C';
        cov_x = A_power_target_time * initial_state.parameters.covariance *...
                A_power_target_time' + flipped_ctrb_mat_disturb *...
                                 cov_concat_disturb * flipped_ctrb_mat_disturb';
    else
        % Computation of mean and covariance of x at target_time by Proposition
        % 1, HSCC 2017
        % mu_x = A^\tau * x_0 + C * mu_W;
        mean_x = A_power_target_time * initial_state + ...
                                flipped_ctrb_mat_disturb * mean_concat_disturb;
        % cov_x = C * cov_{W} * C';
        cov_x = flipped_ctrb_mat_disturb * cov_concat_disturb *...
                                                      flipped_ctrb_mat_disturb';
    end
end
