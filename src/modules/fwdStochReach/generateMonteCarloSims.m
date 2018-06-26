function concat_state_realization = generateMonteCarloSims(...
                                               n_monte_carlo_sims, ...
                                               sys, ...
                                               initial_state, ...
                                               time_horizon, ...
                                               varargin)
% SReachTools/stochasticReachAvoid/generateMonteCarloSims: Generate Monte-Carlo
% simulations for a Gaussian LTI system (controlled or uncontrolled)
% ============================================================================
% generateMonteCarloSims produces a required number of trajectories for a
% Gaussian LTI system.
%
% Usage: See checkViaMonteCarloSims and examples/forwardStochasticReachCWH.mlx
%
% =============================================================================
% concat_state_realization = generateMonteCarloSims(n_monte_carlo_sims, ...
%                                                   sys, ...
%                                                   initial_state, ...
%                                                   time_horizon, ...
%                                                   optimal_input_vector)
%
% Inputs:
% -------
%   n_monte_carlo_sims   - Number of Monte-Carlo simulation particles to be used 
%                          for estimation of the reach-avoid probability
%   sys                  - LtiSystem object describing the system to be verified
%   initial_state        - Deterministic x_0
%   time_horizon         - Time horizon (N) of the stochastic reach-avoid
%                          problem
%   optimal_input_vector - (Optional) Optimal open-loop policy. Required only if 
%                          the system is controlled
%
%
% Outputs:
% --------
%   concat_state_realization- Matrix of concatenate state (row) vectors
%                                   stacked columnwise. Each row comprises of
%                                   the state trajectory as [x_1; x_2; ...; x_N]
%
% See also checkViaMonteCarloSims
%
% Notes:
% ------
% * MATLAB DEPENDENCY: Uses MATLAB's Statistics and Machine Learning Toolbox
%   (mvnrnd)
% * INPUT HANDLING: Delegates part of input handling to @LtiSystem/getConcatMats
% * Assumes IID Gaussian disturbance for the LTI system. 
% * For uncontrolled system, the optimal_input_vector NEED NOT be provided
% * For controlled system, an open-loop controller NEEDS to be provided. The
%   optimal_input_vector should be a ((sys.input_dim) *
%   time_horizon)-dimensional vector U = [u_0; u_1; ...; u_N] (column vector).
% 
% ============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/abyvinod/SReachTools/blob/master/LICENSE
% 
%

    % EXTERNAL DEPENDENCY CHECK
    assert(exist('mvnrnd','file')==2, ...
           'SReachTools:setup_error', ...
           ['This function needs MATLAB''s Statistics and Machine Learning', ...
            ' Toolbox.']);
    %% Input handling 
    % Ensure that n_monte_carlo_sims is a scalar and positive
    assert(isscalar(n_monte_carlo_sims) &&...
           n_monte_carlo_sims > 0, ...
           'SReachTools:invalidArgs', ...
           'Expected a scalar positive n_monte_carlo_sims');
    % Ensure initial state vector is valid or it is a random vector
    if ~isa(initial_state, 'RandomVector')
        assert(isvector(initial_state) &&...
               length(initial_state) == sys.state_dim, ...
               'SReachTools:invalidArgs', ...
               'Expected a valid dimensioned initial_state');
    end
    % Ensure sys has a Gaussian disturbance 
    % GUARANTEES: well-defined sys.dist_dim > 0 and
    %                          sys.dist.parameters.{mean,covariance}
    assert(strcmp(sys.dist.type,'Gaussian'), ...
           'SReachTools:invalidArgs', ...
           ['Monte-Carlo simulations currently only for Gaussian-perturbed', ...
            ' systems']);
    % Compute the concatenated matrices for X
    % GUARANTEES: Scalar time_horizon>0
    [Abar, H, G_matrix] = getConcatMats(sys, time_horizon);
    % Input handling of optimal_input_vector and creation of
    % concat_optimal_input_vector
    if sys.input_dim > 0
        optimal_input_vector = varargin{1};
        % Ensure optimal_input_vector is a valid open loop controller
        assert(isvector(optimal_input_vector) &&...
               length(optimal_input_vector) == ...
                    time_horizon * sys.input_dim, ...
               'SReachTools:invalidArgs', ...
               'Expected a valid dimensioned optimal_input_vector');
        % Repeat the same open-loop policy across all the simulations
        concat_optimal_input_vector = ...
                  repmat(optimal_input_vector,1, n_monte_carlo_sims);
    else
        concat_optimal_input_vector = ...
                zeros(0, n_monte_carlo_sims);
    end

    if isa(initial_state,'RandomVector')
        % Draw from the initial state PDF
        concat_initial_state = mvnrnd(initial_state.parameters.mean', ...
                                      initial_state.parameters.covariance, ...
                                      n_monte_carlo_sims)';
    else
        % Concatenation of initial state vector
        concat_initial_state = repmat(initial_state,1, n_monte_carlo_sims);
    end
    % Compute the concatenated disturbance mean (repeats mean which is a row
    % vectors time_horizon times row-wise)
    concat_disturb_mean = repmat(sys.dist.parameters.mean, ...
                                           time_horizon, ...
                                           1)';
    % Compute the concatenated disturbance covariance (creates a diagonal matrix
    % under the IID assumption)
    concat_disturb_cov = kron(eye(time_horizon), ...
                              sys.dist.parameters.covariance);
    % Realization of the concatenated disturbance random vectors with each
    % realization stored columnwise
    concat_disturb_realizations = mvnrnd(concat_disturb_mean, ...
                                         concat_disturb_cov, ...
                                         n_monte_carlo_sims)';
    % Realization of the random vector X (columnwise)
    % See @LtiSystem/getConcatMats for more info on the notation used
    concat_state_realization = Abar * concat_initial_state +...
                          H * concat_optimal_input_vector +...
                          G_matrix * concat_disturb_realizations;
end
