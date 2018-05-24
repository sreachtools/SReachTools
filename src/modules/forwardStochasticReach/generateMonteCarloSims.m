function concatenated_state_realization = generateMonteCarloSims(...
                                               no_of_monte_carlo_simulations,...
                                               sys,...
                                               initial_state,...
                                               time_horizon,...
                                               varargin)
% SReachTools/stochasticReachAvoid/generateMonteCarloSims: Generate Monte-Carlo
% simulations for a Gaussian LTI system (controlled or uncontrolled)
% ============================================================================
% generateMonteCarloSims produces a required number of trajectories for a
% Gaussian LTI system.
%
% USAGE: See SReachTools/stochasticReachAvoid/checkViaMonteCarloSims for usage
%
% =============================================================================
% concatenated_state_realization = generateMonteCarloSims(
%                                              no_of_monte_carlo_simulations,...
%                                              sys,...
%                                              initial_state,...
%                                              time_horizon,...
%                                              optimal_input_vector)
%
% Inputs:
% -------
%   no_of_monte_carlo_simulations - Number of Monte-Carlo simulation particles
%                                   to be used for estimation of the reach-avoid
%                                   probability
%   sys                           - LtiSystem object describing the system to be
%                                   verified
%   initial_state                 - Initial state of interest
%   time_horizon                  - Time horizon (N) of the stochastic reach-avoid
%                                   problem
%   optimal_input_vector          - (Optional) Optimal open-loop policy.
%                                   Required only if the system is controlled
%
%
% Outputs:
% --------
%   concatenated_state_realization- Matrix of concatenate state (row) vectors
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
%   optimal_input_vector should be a ((sys.input_dimension) *
%   time_horizon)-dimensional vector U = [u_0; u_1; ...; u_N] (column vector).
% 
% ============================================================================
% 
% This function is part of the Stochastic Optimal Control Toolbox.
% License for the use of this function is given in
%      https://github.com/abyvinod/SReachTools/blob/master/LICENSE
% 
%

    % EXTERNAL DEPENDENCY CHECK
    assert(exist('mvnrnd','file')==2,...
           'SReachTools:setup_error',...
           ['This function needs MATLAB''s Statistics and Machine Learning',...
            ' Toolbox.']);
    %% Input handling 
    % Ensure that no_of_monte_carlo_simulations is a scalar and positive
    assert(isscalar(no_of_monte_carlo_simulations) &&...
           no_of_monte_carlo_simulations > 0,...
           'SReachTools:invalidArgs',...
           'Expected a scalar positive no_of_monte_carlo_simulations');
    % Ensure initial state vector is valid
    assert(isvector(initial_state) &&...
           length(initial_state) == sys.state_dimension,...
           'SReachTools:invalidArgs',...
           'Expected a valid dimensioned initial_state');
    % Ensure sys has a Gaussian disturbance 
    % GUARANTEES: well-defined sys.disturbance_dimension > 0 and
    %                          sys.disturbance.parameters.{mean,covariance}
    assert(strcmp(sys.disturbance.type,'Gaussian'),...
           'SReachTools:invalidArgs',...
           ['Monte-Carlo simulations currently only for Gaussian-perturbed',...
            ' systems']);
    % Compute the concatenated matrices for X
    % GUARANTEES: Scalar time_horizon>0
    [Abar, H_matrix, G_matrix] = getConcatMats(sys, time_horizon);
    % Input handling of optimal_input_vector and creation of
    % concatenated_optimal_input_vector
    if sys.input_dimension > 0
        optimal_input_vector = varargin{1};
        % Ensure optimal_input_vector is a valid open loop controller
        assert(isvector(optimal_input_vector) &&...
               length(optimal_input_vector) ==...
                                          time_horizon * sys.input_dimension,...
               'SReachTools:invalidArgs',...
               'Expected a valid dimensioned optimal_input_vector');
        % Repeat the same open-loop policy across all the simulations
        concatenated_optimal_input_vector = ...
                  repmat(optimal_input_vector,1, no_of_monte_carlo_simulations);
    else
        concatenated_optimal_input_vector =...
                zeros(1, no_of_monte_carlo_simulations);
    end

    % Concatenation of initial state vector
    concatenated_initial_state = ...
                         repmat(initial_state,1, no_of_monte_carlo_simulations);
    % Compute the concatenated disturbance mean (repeats mean which is a row
    % vectors time_horizon times row-wise)
    concatenated_disturbance_mean = repmat(sys.disturbance.parameters.mean,...
                                           time_horizon,...
                                           1)';
    % Compute the concatenated disturbance covariance (creates a diagonal matrix
    % under the IID assumption)
    concatenated_disturbance_covariance =...
                                    kron(eye(time_horizon),... 
                                         sys.disturbance.parameters.covariance);
    % Realization of the concatenated disturbance random vectors with each
    % realization stored columnwise
    concatenated_disturbance_realizations =...
                                  mvnrnd(concatenated_disturbance_mean,...
                                         concatenated_disturbance_covariance,...
                                         no_of_monte_carlo_simulations)';
    % Realization of the random vector X (columnwise)
    % See @LtiSystem/getConcatMats for more info on the notation used
    concatenated_state_realization = Abar * concatenated_initial_state +...
                          H_matrix * concatenated_optimal_input_vector +...
                          G_matrix * concatenated_disturbance_realizations;
end
