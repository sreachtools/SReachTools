function [varargout] = getHmatMeanCovForXSansInput(sys, ...
                                                   initial_state, ...
                                                   time_horizon)
% Get input policy-free mean and covariance of the trajectory from a given 
% initial state for a known time horizon and the concatenated input matrix
% ============================================================================
%
% Helps in the computation of the mean and covariance of the concatenated
% state vector X for a given stochastic LTI system as given in (17) of
%
% A. Vinod and M. Oishi, "Scalable Underapproximation for Stochastic Reach-Avoid
% Problem for High-Dimensional LTI Systems using Fourier Transforms," in IEEE
% Control Systems Letters (L-CSS), 2017.
%
% Also, returns H, and Z and G if needed
%
% For more details on the matrix notation, please see the documentation of
% LtvSystem/getConcatMats(). 
%
% Usage: See getLowerBoundStochReachAvoid
%
% ============================================================================
% 
% [H, mean_X_sans_input, cov_X_sans_input, varargout] = ...
%                getHmatMeanCovForXSansInput(sys, ...
%                                            initial_state, ...
%                                            time_horizon)
% Inputs:
% -------
%   sys           - An object of LtvSystem class 
%   initial_state - Initial state can be a deterministic n-dimensional vector
%                   x_0 or a RandomVector object
%   time_horizon  - Time of interest (N)
%
% Outputs:
% --------
%   H                - Concatenated input matrix
%   mean_X_sans_input- Mean of X with zero input under the disturbance from the
%                      provided initial state
%   cov_X_sans_input - Covariance of X with zero input under the disturbance
%                      from the provided initial state
%   Z                - (optional) Concatenated state matrix
%   G                - (optional) Concatenated disturbance matrix
%
% Notes:
% ------
% * X refers to the concatenated state vector X=[x_1^\top x_2^\top ...
%   x_N^\top]^\top. See @LtvSystem/getConcatMats for more
%   information about the notation used.
% * Assumes the disturbance is independent and identically distributed
% * This function also serves as a delegatee for input handling.
% 
% ============================================================================
%
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
%
%

    %% Input handling
    inpar = inputParser();
    inpar.addRequired('sys', @(x) validateattributes(x,...
        {'LtiSystem','LtvSystem'}, {'nonempty'}));
    % Ensure that the initial_state is a 
    inpar.addRequired('initial_state', @(x) validateattributes(x, ...
        {'RandomVector', 'numeric'}, {'nonempty'}))
    % Ensure that the target_time is a positive scalar
    inpar.addRequired('time_horizon', @(x) validateattributes(x, ...
        {'numeric'}, {'scalar', 'integer', '>', 0}));

    try
        inpar.parse(sys, initial_state, time_horizon);
    catch err
        exc = SrtInvalidArgsError.withFunctionName();
        exc = exc.addCause(err);
        throwAsCaller(exc);
    end

    % Ensure that initial state is a column vector of appropriate dimension OR
    % random vector of approriate dimension
    if isa(initial_state,'RandomVector') &&...
        (initial_state.dim ~= sys.state_dim ||...
         ~strcmpi(initial_state.type, 'Gaussian'))
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

    %% Compute the concatenated matrices for X
    % H will be a zeros(sys.state_dim * time_horizon, 1) matrix for an
    % uncontrolled LTI system
    [Z, H, G] = getConcatMats(sys, time_horizon);

    % IID assumption allows to compute the mean and covariance of the
    % concatenated disturbance vector W
    mean_concat_disturb = kron(ones(time_horizon,1), ...
                               sys.dist.parameters.mean);
    cov_concat_disturb  = kron(eye(time_horizon), ...
                               sys.dist.parameters.covariance);
                                 

    if isa(initial_state,'RandomVector')
        % TODO: Waiting for a reference --- but essentially propagation of mean
        % and covariance
        mean_X_sans_input = Z * initial_state.parameters.mean + G *...
            mean_concat_disturb;
        cov_X_sans_input = Z * initial_state.parameters.covariance * Z' + ...
            G * cov_concat_disturb * G';
    else
        % Computation of mean and covariance of X (sans input) by (17),LCSS 2017
        mean_X_sans_input = Z * initial_state + G * mean_concat_disturb;
        cov_X_sans_input = G * cov_concat_disturb * G';
    end
    
    % Ensure that cov_X_sans_input is a symmetric matrix
    if ~issymmetric(cov_X_sans_input)
        % Compute the symmetric component of it
        symm_cov_X = (cov_X_sans_input+cov_X_sans_input')/2;
        % Max error element-wise
        max_err = max(max(abs(cov_X_sans_input - symm_cov_X)));
        if max_err > eps
            warning(sprintf(['Non-symmetric covariance matrix made ',...
                'symmetric (max elementwise error: %1.3e)!'], max_err));
        end
        cov_X_sans_input = symm_cov_X;
    end
    
    % Optional output arguments
    varargout{1} = H;
    varargout{2} = mean_X_sans_input;
    varargout{3} = cov_X_sans_input;
    varargout{4} = Z;
    varargout{5} = G;
end
