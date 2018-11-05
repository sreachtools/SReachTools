function [concat_state_realization, ...
          concat_disturb_realizations]= generateMonteCarloSims(...
                                               n_monte_carlo_sims, ...
                                               sys, ...
                                               initial_state, ...
                                               time_horizon, ...
                                               varargin)
% Generate Monte-Carlo simulations for a Gaussian-perturbed LTI/LTV system
% (controlled or uncontrolled)
% ============================================================================
% 
% generateMonteCarloSims produces a required number of trajectories for a
% Gaussian LTI system.
%
% See also examples/forwardStochasticReachCWH.m
%
% =============================================================================
% [concat_state_realization, concat_disturb_realizations] = ...
%    generateMonteCarloSims(n_monte_carlo_sims, sys, initial_state, ...
%       time_horizon, optimal_input_vector, optimal_input_gain)
%
% Inputs:
% -------
%   n_monte_carlo_sims   - Number of Monte-Carlo simulation particles to be used 
%                          for estimation of the reach-avoid probability
%   sys                  - System description as a LtiSystem/LtvSystem object
%   initial_state        - Deterministic x_0
%   time_horizon         - Time horizon (N) of the stochastic reach-avoid
%                          problem
%   optimal_input_vector - [Optional] Open-loop controller, a column vector of
%                          dimension (sys.input_dim*N) x 1 | Required only if 
%                          the system is controlled
%   optimal_input_gain   - [Optional] Affine disturbance feedback gain, a matrix
%                          of dimension (sys.input_dim*N) x (sys.dist_dim*N) |
%                          Required only if the system is controlled and the
%                          controller is affine disturbance feedback
%
% Outputs:
% --------
%   concat_state_realization  - Matrix of concatenate state (row) vectors
%                               stacked columnwise. Each row comprises of the
%                               state trajectory as [x_0; x_1; x_2; ...; x_N]
%   concat_disturb_realization- Matrix of concatenate disturbance (row) vectors
%                               stacked columnwise. Each row comprises of
%                               the state trajectory as [w_0; w_1; ...; w_{N-1}]
%
% Notes:
% ------
% * Assumes IID Gaussian disturbance for the LTI/LTV system. 
% * For controlled system, an open-loop controller NEEDS to be provided. The
%   optimal_input_vector should be a ((sys.input_dim) * time_horizon)-dim.
%   vector U = [u_0; u_1; ...; u_N] (column vector).
% * For uncontrolled system, the optimal_input_vector NEED NOT be provided
% 
% ============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 
%

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
    [Z, H, G] = getConcatMats(sys, time_horizon);
    % Input handling of optimal_input_vector and creation of
    % concat_optimal_input_vector
    concat_optimal_input_vector = zeros(time_horizon * sys.input_dim,...
        n_monte_carlo_sims);
    optimal_input_dist_gain = zeros(time_horizon * sys.input_dim,...
        time_horizon * sys.dist_dim);
    if sys.input_dim > 0
        if nargin >= 5
            optimal_input_vector = varargin{1};
            % Ensure optimal_input_vector is a valid open loop controller
            assert(isvector(optimal_input_vector) &&...
                   length(optimal_input_vector)==time_horizon*sys.input_dim, ...
                   'SReachTools:invalidArgs', ...
                   'Expected a valid dimensioned optimal_input_vector');
            % Repeat the same open-loop policy across all the simulations
            concat_optimal_input_vector = ...
                      repmat(optimal_input_vector,1, n_monte_carlo_sims);
            if nargin == 6
                optimal_input_dist_gain = varargin{2};
                % Ensure optimal_input_vector is a valid open loop controller
                assert(ismatrix(optimal_input_dist_gain) &&...
                       all(size(optimal_input_dist_gain) == ...
                            [time_horizon * sys.input_dim,...
                             time_horizon * sys.dist_dim]), ...
                       'SReachTools:invalidArgs', ...
                       'Expected a valid dimensioned optimal_input_dist_gain');                
            end
        end
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
    W = sys.dist.concat(time_horizon);
    % Realization of the concatenated disturbance random vectors with each
    % realization stored columnwise
    concat_disturb_realizations = mvnrnd(W.parameters.mean, ...
                                         W.parameters.covariance, ...
                                         n_monte_carlo_sims)';
    % Realization of the random vector X (columnwise)
    % See @LtiSystem/getConcatMats for more info on the notation used
    concat_state_realization = [concat_initial_state;
        Z * concat_initial_state + H * concat_optimal_input_vector + ...
          (H * optimal_input_dist_gain + G) * concat_disturb_realizations];
end
