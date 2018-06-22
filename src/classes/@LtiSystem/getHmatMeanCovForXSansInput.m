function [H, mean_X_sans_input, cov_X_sans_input, varargout] = ...
                            getHmatMeanCovForXSansInput(sys, ...
                                                        initial_state, ...
                                                        time_horizon)
% SReachTools/LtiSystem/getHmatMeanCovForXSansInput: Get input policy-free mean
% and covariance of the trajectory from a given initial state for a known time
% horizon and the concatenated input matrix
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
% LtiSystem/getConcatMats(). 
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
%   sys           - An object of LtiSystem class 
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
%   x_N^\top]^\top. See @LtiSystem/getConcatMats for more
%   information about the notation used.
% * This function also serves as a delegatee for input handling.
% 
% ============================================================================
%
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/abyvinod/SReachTools/blob/master/LICENSE
%
%

    %% Input handling
    % Ensure that the given system has a Gaussian disturbance
    assert( (strcmp(class(sys.dist),'StochasticDisturbance') ||...
             strcmp(class(sys.dist),'RandomVector')) &&...
            strcmp(sys.dist.type,'Gaussian'), ...
           'SReachTools:invalidArgs', ...
           ['getHmatMeanCovForXSansInput is for', ...
            ' Gaussian-perturbed LTI systems only']);

    %% Compute the concatenated matrices for X
    % GUARANTEES: Scalar time_horizon>0
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
        mean_X_sans_input = Z * initial_state.parameters.mean + G *...
            mean_concat_disturb;
        cov_X_sans_input = Z * initial_state.parameters.covariance * Z' + ...
            G * cov_concat_disturb * G';
    else
        % Ensure that initial state is a column vector of appropriate dimension
        assert( size(initial_state,1) == sys.state_dim &&...
                size(initial_state,2) == 1, ...
               'SReachTools:invalidArgs', ...
               ['Expected a sys.state_dim-dimensional column-vector ', ...
                'for initial state']);

        % Computation of mean and covariance of X (sans input) by (17), LCSS 2017
        mean_X_sans_input = Z * initial_state + G * mean_concat_disturb;
        cov_X_sans_input = G * cov_concat_disturb * G';
    end

    % Optional output arguments
    varargout{1} = Z;
    varargout{2} = G;
end
