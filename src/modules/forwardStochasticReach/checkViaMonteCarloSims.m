function [reach_avoid_probability_mcarlo,...
          legend_cell] = checkViaMonteCarloSims(...
                                               no_of_monte_carlo_simulations,...
                                               sys,...
                                               initial_state,...
                                               time_horizon,...
                                               safe_set,...
                                               target_set,...
                                               optimal_input_vector,...
                                               legend_cell,...
                                               no_of_simulations_to_plot)
% SReachTools/forwardStochasticReach/checkViaMonteCarloSims: Monte-Carlo
% simulations-based probability computation of lying in a target tube
% ============================================================================
% checkViaMonteCarloSims uses Monte-Carlo simulations to approximate the
% terminal hitting-time stochastic reach-avoid problem from an initial
% state given an optimal open-loop input policy.
%
% USAGE: See examples/verificationCwhForAnInitialState.m for usage
%
% =============================================================================
% [reach_avoid_probability_mcarlo,...
%  legend_cell] = checkViaMonteCarloSims(...
%                                              no_of_monte_carlo_simulations,...
%                                              sys,...
%                                              initial_state,...
%                                              time_horizon,...
%                                              safe_set,...
%                                              target_set,...
%                                              optimal_input_vector,...
%                                              legend_cell,...
%                                              no_of_simulations_to_plot)
% 
% Inputs:
% -------
%   no_of_monte_carlo_simulations- Number of Monte-Carlo simulation particles
%                                  to be used for estimation of the reach-avoid
%                                  probability
%   sys                          - LtiSystem object describing the system to be
%                                  verified
%   initial_state                - Initial state of interest
%   time_horizon                 - Time horizon of the stochastic reach-avoid
%                                  problem
%   safe_set                     - Safe set for stochastic reach-avoid problem
%   target_set                   - Target set for stochastic reach-avoid problem
%   optimal_input_vector         - Optimal open-loop policy
%                                  ((sys.input_dimension) *
%                                  time_horizon)-dimensional vector 
%                                  U = [u_0; u_1; ...; u_N] (column vector)
%   legend_cell                  - Cell array containing legend strings (will
%                                  be updated by this function if plotting is
%                                  enabled)
%   no_of_simulations_to_plot    - Plot how many simulations from this initial
%                                  state
%
% Outputs:
% --------
%   reach_avoid_probability_mcarlo- Terminal-hitting stochastic reach avoid
%                                   probability via Monte-Carlo simulation
%   legend_cell                   - Cell array containing legend strings (will
%                                   be updated by this function if plotting is
%                                   enabled)
%
% See also generateMonteCarloSims,
% getFtLowerBoundStochasticReachAvoid,
% @LtiSystem/getConcatMats
%
% Notes:
% ------
% * NOT ACTIVELY TESTED: Builds on other tested functions.
% * MATLAB DEPENDENCY: Uses MATLAB's Statistics and Machine Learning Toolbox
%   (mvnrnd)
% * INPUT HANDLING: Delegates part of input handling to
%   generateMonteCarloSims
% * Assumes IID Gaussian disturbance for the LTI system. 
% * For uncontrolled system, the optimal_input_vector needs to be empty.
% * For controlled system, an open-loop controller needs to be provided.
% * See @LtiSystem/getConcatMats for more information about the
%   notation used.
% 
% ============================================================================
% 
% This function is part of the Stochastic Optimal Control Toolbox.
% License for the use of this function is given in
%      https://github.com/abyvinod/SReachTools/blob/master/LICENSE
% 
%

    %% Input handling
    % Construct concatenated target tube (Called earlier to ensure
    % safe_set.contains will work)
    % GUARANTEES: Non-empty target and safe sets (polyhedron) and scalar
    %             time_horizon>0
    [concatenated_target_tube_A, concatenated_target_tube_b] =...
                                        getConcatTargetTube(safe_set,...
                                                                  target_set,...
                                                                  time_horizon);
    % Ensure that legend_cell is a cell array
    assert(iscell(legend_cell),...
           'SReachTools:invalidArgs',...
           'Expected a cell array for legend_cell');
    % Ensure that time horizon is a scalar and positive
    assert(isscalar(no_of_simulations_to_plot) &&...
           no_of_simulations_to_plot >= 0,...
           'SReachTools:invalidArgs',...
           'Expected a scalar non-negative no_of_simulations_to_plot');

    % Get no_of_monte_carlo_simulations number of trajectories under the
    % provided open-loop policy and initial state
    % GUARANTEES: well-defined no_of_monte_carlo_simulations, initial_state,
    %                          time_horizon, optimal_input_vector
    concatenated_state_realization = generateMonteCarloSims(...
                                               no_of_monte_carlo_simulations,...
                                               sys,...
                                               initial_state,...
                                               time_horizon,...
                                               optimal_input_vector);

    % Define a boolean vector that is used for relative frequency estimation of
    % the reach-avoid probability
    pass_fail_vector = zeros(no_of_monte_carlo_simulations,1);
    % For plotting
    if no_of_simulations_to_plot > 0
        plot_at_index_intervals =...
                 round(no_of_monte_carlo_simulations/no_of_simulations_to_plot);
        green_legend_updated = 0;
        red_legend_updated = 0;
    end
    % For MonteCarlo simulation and plotting and realization-by-realization
    % check of satisfaction of the reach-avoid objective
    for realization_index = 1: no_of_monte_carlo_simulations
        % Check if the trajectory satisfies the reach-avoid objective
        if concatenated_target_tube_A * ...
                       concatenated_state_realization(:,realization_index) <=...
                                concatenated_target_tube_b
            % realization is reach-avoid satisfying
            pass_fail_vector(realization_index) = 1;
            % Assign green triangle as the marker
            markerString = 'g^-';
        else
            % Assign red asterisk as the marker
            markerString = 'r*-';
        end
        % Plotting
        if no_of_simulations_to_plot > 0 
            if abs(mod(realization_index, plot_at_index_intervals)) < 1e-8
                % Create [x(t_1) x(t_2)... x(t_N)]
                reshaped_X_vector = ...
                 reshape(concatenated_state_realization(:,realization_index),...
                        sys.state_dimension,[]);
                % This realization is to be plotted
                h = plot([initial_state(1), reshaped_X_vector(1,:)],...
                         [initial_state(2), reshaped_X_vector(2,:)],...
                         markerString, 'MarkerSize',10);
                % Update the legends if the first else, disable
                if strcmp(markerString,'g^-')
                    if green_legend_updated
                        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    else
                        green_legend_updated = 1;
                        legend_cell{end+1} = 'Good trajectory';
                    end
                elseif strcmp(markerString,'r*-')
                    if red_legend_updated
                        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    else
                        red_legend_updated = 1;
                        legend_cell{end+1} = 'Bad trajectory';
                    end
                end
            end
        end
    end
    % Compute the reach_avoid_probability_mcarlo via relative frequency
    reach_avoid_probability_mcarlo =...
                sum(pass_fail_vector)/no_of_monte_carlo_simulations;
end
