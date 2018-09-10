function [mean_x, cov_x] = getFSRPDMeanCov(sys, ...
                                           initial_state, ...
                                           target_time)
% Compute the mean and the covariance of the state at a time instant in future
% ============================================================================
%
% Computes the mean and the covariance of a Gaussian-perturbed LTI/LTV
% uncontrolled system. This function implements Proposition 1 of
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
%   sys           - An object of LtiSystem/LtvSystem class 
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
    % input handling overlapping with getHmatMeanCovForXSansInput with
    % time_horizon replaced with target_time (positivity relaxed)
    inpar = inputParser();
    inpar.addRequired('sys', @(x) validateattributes(x,...
        {'LtiSystem','LtvSystem'}, {'nonempty'}));
    % Ensure that the initial_state is a 
    inpar.addRequired('initial_state', @(x) validateattributes(x, ...
        {'RandomVector', 'numeric'}, {'nonempty'}))
    % Ensure that the target_time is a positive scalar
    inpar.addRequired('target_time', @(x) validateattributes(x, ...
        {'numeric'}, {'scalar', 'integer', '>=', 0}));

    try
        inpar.parse(sys, initial_state, target_time);
    catch err
        exc = SrtInvalidArgsError.withFunctionName();
        exc = exc.addCause(err);
        throwAsCaller(exc);
    end

    % Ensure that initial state is a column vector of appropriate dimension OR
    % random vector of approriate dimension
    if isa(initial_state,'RandomVector') &&...
        initial_state.dim ~= sys.state_dim &&...
        strcmp(initial_state.type, 'Gaussian')
        exc = SrtInvalidArgsError.withFunctionName();
        exc = exc.addCause(SrtInvalidArgsError(['Expected a ', ...
            'sys.state_dim-dimensional Gaussian random vector for initial ',...
            'state']));
        throw(exc);
    elseif isa(initial_state,'numeric') &&...
        ~isequal(size(initial_state), [sys.state_dim 1])
        exc = SrtInvalidArgsError.withFunctionName();
        exc = exc.addCause(SrtInvalidArgsError(['Expected a ', ...
            'sys.state_dim-dimensional column-vector for initial state']));
        throw(exc);
    end

    % Ensure that the given system has a Gaussian disturbance
    if isa(sys.dist, 'StochasticDisturbance') || isa(sys.dist, 'RandomVector')
        if ~strcmpi(sys.dist.type, 'Gaussian')
            exc = SrtInvalidArgsError.withFunctionName();
            exc = exc.addCause(SrtInvalidArgsError(['Expected a ',...
                'Gaussian-perturbed LTI system']));
            throw(exc);
        end
    else
        exc = SrtInvalidArgsError.withFunctionName();
        exc = exc.addCause(SrtInvalidArgsError(['Expected a ',...
            'stochastic LTI system']));
        throw(exc);
    end

    %% inputHandling specific to getFSRPDMeanCov
    % Ensure that the given system is uncontrolled
    if sys.input_dim ~= 0
        exc = SrtInvalidArgsError.withFunctionName();
        exc = exc.addCause(SrtInvalidArgsError(['getFSRPDMeanCov ', ...
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

    %% Obtained from getHmatMeanCovForXSansInput
    % Ensure that cov_x is a symmetric matrix
    if ~issymmetric(cov_x)
        % Compute the symmetric component of it
        symm_cov_x = (cov_x+cov_x')/2;
        % Max error element-wise
        max_err = max(max(abs(cov_x - symm_cov_x)));
        if max_err > eps
            warning(sprintf(['Non-symmetric covariance matrix made ',...
                'symmetric (max elementwise error: %1.3e)!'], max_err));
        end
        cov_x = symm_cov_x;
    end
end
