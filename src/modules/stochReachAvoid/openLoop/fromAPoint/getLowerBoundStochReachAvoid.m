function [lb_stoch_reach_avoid, optimal_input_vector] =...
                           getLowerBoundStochReachAvoid(sys,...
                                                        initial_state,...
                                                        target_tube,...
                                                        init_safe_set,...
                                                        varargin)
% SReachTools/stochasticReachAvoid/getLowerBoundStochReachAvoid: Solve the
% stochastic reach-avoid problem (lower bound on the probability and an
% open-loop controller synthesis) using Fourier transform and convex
% optimization
% =============================================================================
%
% getLowerBoundStochReachAvoid implements the Fourier transform-based
% underapproximation to the terminal hitting-time stochastic reach-avoid problem
% discussed in
%
% A. Vinod and M. Oishi, "Scalable Underapproximation for Stochastic Reach-Avoid
% Problem for High-Dimensional LTI Systems using Fourier Transforms," in IEEE
% Control Systems Letters (L-CSS), 2017.
%
% This function has input handling and sets up the problem to be solved by
% computeFtLowerBoundStochReachAvoid.
%
% USAGE: See examples/verificationOfCwhDynamicsForAnInitialState.m
%
% =============================================================================
%
% [lb_stoch_reach_avoid, optimal_input_vector] =...
%                                 getLowerBoundStochReachAvoid(sys,...
%                                                              initial_state,...
%                                                              target_tube,...
%                                                              varargin)
%
% Inputs:
% -------
%   sys                  - LtiSystem object describing the system to be verified
%   initial_state        - Initial state of interest
%   target_tube          - Target tube to stay within [TargetTube object]
%   init_safe_set        - Safe set for initial state
%   guess_optimal_input_vector
%                        - (Optional) Provide a concatenated guess for the
%                          optimal input policy vector in the form of U = [u_0;
%                          u_1; ...; u_N]. [If unsure, provide []. This will
%                          trigger a CVX-based initialization computation.]
%   desired_accuracy     - (Optional) Accuracy  [Default 5e-3]
%   PSoptions            - (Optional) Options for patternsearch [Default
%                           psoptimset('Display', 'off')]
%
% Outputs:
% --------
%   lb_stoch_reach_avoid - Lower bound on the terminal-hitting stochastic reach
%                          avoid problem computed using Fourier transform and
%                          convex optimization
%   optimal_input_vector - Optimal open-loop policy ((sys.input_dimension) *
%                          time_horizon)-dim.  vector U = [u_0; u_1; ...; u_N]
%                          (column vector)
%
% See also computeFtLowerBoundStochReachAvoid, getCcLowerBoundStochReachAvoid.
%
% Notes:
% * NOT ACTIVELY TESTED: Builds on other tested functions.
% * MATLAB DEPENDENCY: Uses MATLAB's Global Optimization Toolbox; Statistics and
%                      Machine Learning Toolbox.
%                      Needs patternsearch for gradient-free optimization
%                      Needs normpdf, normcdf, norminv for Genz's algorithm
% * EXTERNAL DEPENDENCY: Uses MPT3 and CVX (optional)
%                        Needs MPT3 for defining a controlled system and the
%                        definition of the safe and the target (polytopic) sets
%                        Needs CVX to setup a convex optimization problem that
%                        initializes the patternsearch-based optimization. If 
%                        CVX is unavailable, the user may provide a guess for 
%                        the initialization.
% * Specify both desired_accuracy and PSoptions or neither to use the defaults 
% * Specify an optional guess_optimal_input_vector to skip the use of CVX
% * See @LtiSystem/getConcatMats for more information about the
%     notation used.
% 
% =============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/abyvinod/SReachTools/blob/master/LICENSE
%
%

    % Get half space representation of the target tube and time horizon
    [concat_target_tube_A, concat_target_tube_b] = target_tube.concat();
    time_horizon = length(target_tube);

    % Input handling: Check for safe set
    validateattributes(init_safe_set, {'Polyhedron'}, {'nonempty'});


    % Check if safe set contains the initial state
    if ~init_safe_set.contains(initial_state)
        % Stochastic reach-avoid probability is zero and no admissible open-loop
        % policy exists, if given an unsafe initial state
        lb_stoch_reach_avoid = 0;
        optimal_input_vector = nan(sys.input_dimension * time_horizon, 1);
    else
        % Stochastic reach-avoid probability may be non-trivial
        % EXTERNAL DEPENDENCY CHECK
        assert(exist('cvx_begin','file')==2, ...
               'SReachTools:setup_error', ...
               'This function uses CVX. Please get it from http://cvxr.com.');
        assert(exist('patternsearch','file')==2, ...
               'SReachTools:setup_error', ...
               'This function needs MATLAB''s Global Optimization Toolbox.');
        assert(exist('normcdf','file')==2, ...
               'SReachTools:setup_error', ...
               ['This function needs MATLAB''s Statistics and Machine ', ...
                'Learning Toolbox.']);

        % Construct U^N 
        % GUARANTEES: Non-empty input sets (polyhedron) and scalar
        %             time_horizon>0
        [concat_input_space_A, concat_input_space_b] = ...
                                              getConcatInputSpace(sys, ...
                                                                  time_horizon);
        % Compute H, mean_X_sans_input, cov_X_sans_input for the
        % safety_cost_function definition
        % GUARANTEES: Gaussian-perturbed LTI system (sys) and well-defined
        % initial_state and time_horizon
        [H, mean_X_sans_input, cov_X_sans_input] = ...
                                  getHmatMeanCovForXSansInput(sys, ...
                                                              initial_state, ...
                                                              time_horizon);

        % Parsing the optional arguments 
        if length(varargin) == 3
            % First optional argument is the guess_optimal_input_vector
            if ~isempty(varargin{1})
                assert(isvector(varargin{1}) &&...
                       length(varargin{1}) == ...
                       sys.input_dimension * time_horizon, ...
                       'SReachTools:invalidArgs', ...
                       ['Expected a well-dimensioned row vector ', ...
                        'guess_optimal_input_vector']);
                guess_optimal_input_vector = varargin{1};
            else
                guess_optimal_input_vector = [];
            end
            % Second optional argument is the desired_accuracy
            assert(isscalar(varargin{2}), ...
                   'SReachTools:invalidArgs', ...
                   'Expected a scalar value for desired_accuracy');
            desired_accuracy = varargin{2};
            % Third optional argument is the options for patternsearch,
            % PSoptions (TODO: No validation being done here)
            PSoptions = varargin{3};
        elseif length(varargin) == 1
            % First optional argument is the guess_optimal_input_vector
            if ~isempty(varargin{1})
                assert(isvector(varargin{1}) &&...
                       length(varargin{1}) == ...
                       sys.input_dimension * time_horizon, ...
                       'SReachTools:invalidArgs', ...
                       ['Expected a well-dimensioned row vector ', ...
                        'guess_optimal_input_vector']);
                guess_optimal_input_vector = varargin{1};
            else
                guess_optimal_input_vector = [];
            end
            desired_accuracy = 1e-3;
            PSoptions = psoptimset('Display', 'off');
        elseif isempty(varargin) == 0
            guess_optimal_input_vector = [];
            desired_accuracy = 1e-3;
            PSoptions = psoptimset('Display', 'off');
        else
            error('SReachTools:invalidArgs', ...
                  ['guess_optimal_input_vector, desired_accuracy, and ', ...
                   'PSoptions are the only additional options.']);
        end
        % END OF INPUT HANDLING

        % Patternsearch and Fourier transform-based open-loop underapproximation
        % of the stochastic reach-avoid problem
        [lb_stoch_reach_avoid, optimal_input_vector] = ...
            computeFtLowerBoundStochReachAvoid(sys, ...
                                               initial_state, ...
                                               time_horizon, ...
                                               concat_input_space_A, ... 
                                               concat_input_space_b, ...
                                               concat_target_tube_A, ... 
                                               concat_target_tube_b, ...
                                               H, ...
                                               mean_X_sans_input, ...
                                               cov_X_sans_input, ...
                                               guess_optimal_input_vector, ...
                                               desired_accuracy, ...
                                               PSoptions);
    end
end

%% If an open_loop policy is desired arranged in increasing time columnwise
% optimal_open_loop_control_policy = reshape(optimal_input_vector, ...
%                                    sys.input_dimension, ...
%                                    time_horizon);


%% Patternsearch other options
% mesh_tolerance_for_patternsearch = 1e-6;
% constraint_tolerance_for_patternsearch = 1e-6;
% PSoptions = psoptimset('Display', display_string, ...
%                 'TolMesh', mesh_tolerance_for_patternsearch, ...
%                 'TolCon', constraint_tolerance_for_patternsearch);
