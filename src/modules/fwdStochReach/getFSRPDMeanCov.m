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

    % Input handling
    inpar = inputParser();
    inpar.addRequired('sys', @(x) validateattributes(x, {'LtiSystem'}, ...
        {'nonempty'}));
    inpar.addRequired('initial_state', @(x) validateattributes(x, ...
        {'RandomVector', 'numeric'}, {'nonempty'}))
    
    % Ensure that the target_time is a positive scalar
    inpar.addRequired('target_time', @(x) validateattributes(x, ...
        {'numeric'}, {'scalar', 'integer', '>', 0}));

    try
        inpar.parse(sys, initial_state, target_time);
    catch err
        exc = SrtInvalidArgsError.withFunctionName();
        exc = exc.addCause(err);
        throwAsCaller(exc);
    end

    if ~isa(initial_state,'RandomVector')
        % Ensure that initial state is a column vector of appropriate dimension
        if ~iscolumn(initial_state) && ...
           length(initial_state) == sys.state_dimension

            exc = SrtInvalidArgsError.withFunctionName();
            exc = exc.addCause(SrtInvalidArgsError(['Expected a ', ...
                'sys.state_dimension-dimensional column-vector for ', ...
                'initial state']));
            throw(exc);
        end
    end

    % Ensure that the given system has a Gaussian disturbance
    if ~isa(sys.disturbance), 'StochasticDisturbance') || ...
       ~strcmp(sys.disturbance.type), 'Gaussian')
    
        exc = SrtInvalidArgsError.withFunctionName();
        exc = exc.addCause(SrtInvalidArgsError(['getFSRPDMeanCovariance ', ...
            'accepts only Gaussian-perturbed LTI systems']));
        throw(exc);
    end

    % Ensure that the given system is uncontrolled
    if sys.input_dimension ~= 0
        exc = SrtInvalidArgsError.withFunctionName();
        exc = exc.addCause(SrtInvalidArgsError(['getFSRPDMeanCovariance ', ...
            'accepts only uncontrolled LTI systems']));
        throw(exc);
    end

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
