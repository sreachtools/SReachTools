function [mean_x, cov_x] = getFSRPDMeanCov(sys, ...
                                           initial_state, ...
                                           target_time)
% SReachTools/forwardStochasticReach/getFSRPDMeanCov: Compute the mean and the
% covariance of the state at a time instant in future
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
% Usage: See getProbReachSet, examples/forwardStochasticReachCWH.mlx.
%
% ============================================================================
% 
% [mean_x, cov_x] = getFSRPDMeanCov(sys, initial_state, target_time)
% 
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
% See also getProbReachSet, getHmatMeanCovForXSansInput
%
% Notes:
% ------
% * getHmatMeanCovForXSansInput computes the FSRPD for the joint state vector
%   and getFSRPDMeanCov computes the FSRPD for the state at time t.
%
% ============================================================================
%
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
%
%

    %% Input handling
    if ~isa(initial_state,'RandomVector')
        % Ensure that initial state is a column vector of appropriate dimension
        assert( size(initial_state,1) == sys.state_dim &&...
                size(initial_state,2) == 1, ...
               'SReachTools:invalidArgs', ...
               ['Expected a sys.state_dim-dimensional column-vector ', ...
               'for initial state']);
    end

    % Ensure that the given system has a Gaussian disturbance
    assert( isa(sys.dist,'StochasticDisturbance') ||...
                isa(sys.dist,'RandomVector') &&...
           strcmp(sys.dist.type,'Gaussian') &&...
           sys.input_dim == 0, ...
          'SReachTools:invalidArgs', ...
          'getFSRPDMeanCov accepts only uncontrolled Gaussian-perturbed systems');

    % Ensure that the given system is uncontrolled
    assert(sys.input_dim == 0, ...
          'SReachTools:invalidArgs', ...
          'getFSRPDMeanCov accepts only uncontrolled LTI systems');

    % Ensure that the target_time is a positive scalar
    assert( isscalar(target_time) && target_time > 0, ...
           'SReachTools:invalidArgs', ...
           'Expected a scalar positive target_time');

    % IID assumption allows to compute the mean and covariance of the
    % concatenated disturbance vector W
    mean_concat_disturb = kron(ones(target_time,1), ...
                               sys.dist.parameters.mean);
    cov_concat_disturb =kron(eye(target_time), ...
                             sys.dist.parameters.covariance);

    % Compute the state_transition_matrix and controllability matrix for disturbance
    [Z,~,G] = sys.getConcatMats(target_time);
    state_transition_mat = Z(end-sys.state_dim+1:end,end-sys.state_dim+1:end);
    flipped_ctrb_mat_disturb = G(end-sys.state_dim+1:end,:);

    if isa(initial_state,'RandomVector')
        % mu_x = A^\tau * mu_{x_0} + C * mu_{W};
        mean_x = state_transition_mat * initial_state.parameters.mean + ...
                                flipped_ctrb_mat_disturb * mean_concat_disturb;
        % cov_x = A^\tau * cov_{x_0} * (A^\tau)' + C * cov_{W} * C';
        cov_x = state_transition_mat * initial_state.parameters.covariance *...
                state_transition_mat' + flipped_ctrb_mat_disturb *...
                                 cov_concat_disturb * flipped_ctrb_mat_disturb';
    else
        % Computation of mean and covariance of x at target_time by Proposition
        % 1, HSCC 2017
        % mu_x = A^\tau * x_0 + C * mu_W;
        mean_x = state_transition_mat * initial_state + ...
                                flipped_ctrb_mat_disturb * mean_concat_disturb;
        % cov_x = C * cov_{W} * C';
        cov_x = flipped_ctrb_mat_disturb * cov_concat_disturb *...
                                                      flipped_ctrb_mat_disturb';
    end
end
