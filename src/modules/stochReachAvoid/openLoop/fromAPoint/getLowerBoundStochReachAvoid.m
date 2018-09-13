function [lb_stoch_reach_avoid, optimal_input_vector, varargout] =...
                           getLowerBoundStochReachAvoid(sys,...
                                                        initial_state,...
                                                        target_tube,...
                                                        method,...
                                                        varargin)
% Solve the stochastic reach-avoid problem (lower bound on the probability and 
% an open-loop controller synthesis) using Fourier transform and convex
% optimization
% =============================================================================
%
% getLowerBoundStochReachAvoid computes a lower bound to the terminal
% hitting-time stochastic reach-avoid problem by searching over the space of
% open-loop controllers.
%
% A. Vinod and M. Oishi, "Scalable Underapproximation for Stochastic Reach-Avoid
% Problem for High-Dimensional LTI Systems using Fourier Transforms," in IEEE
% Control Systems Letters (L-CSS), 2017.
%
% This function has input handling and sets up the problem to be solved by
% depending on the user-preference:
%   1. computeFtLowerBoundStochReachAvoid, (Option: genzps)
%   2. computeCcLowerBoundStochReachAvoidIterRisk, (Option: ccciter)
%   3. computeCcLowerBoundStochReachAvoidPwlRisk (Option: cccpwl).
%
% USAGE: TODO
%
% =============================================================================
%
% [lb_stoch_reach_avoid, optimal_input_vector] =...
%                                 getLowerBoundStochReachAvoid(sys,...
%                                                              initial_state,...
%                                                              target_tube,...
%                                                              method,...
%                                                              varargin)
%
% Inputs:
% -------
%   sys              - LtiSystem object
%   initial_state    - Initial state
%   target_tube      - TargetTube object
%   method           - Method to compute the reach-avoid probability
%                         'genzps'        -- Genz's algorithm + Patternsearch
%                         'cccpwl'        -- Piecewise-linear conservative
%                                            implementation for convex
%                                            chance-constrained reformulation
%                         'ccciter'       -- Iterative risk allocation approach 
%                                            for convex chance-constrained
%                                            reformulation
%                         'cccpwl-closed' -- Piecewise-linear conservative
%                                            implementation for convex
%                                            chance-constrained reformulation
%                      See Notes for dependencies.
%   guess_optimal_input_vector
%                    - (Optional) Provide a concatenated guess for the optimal
%                      input policy vector in the form of U = [u_0; u_1; ...;
%                      u_N] for patternsearch. [Default []. This will trigger a
%                      CVX-based computation for initialization of
%                      patternsearch.]
%   desired_accuracy - (Optional) Accuracy for patternsearch  [Default 5e-3]
%   PSoptions        - (Optional) Options for patternsearch [Default
%                       psoptimset('Display', 'off')]
%
% Outputs:
% --------
%   lb_stoch_reach_avoid - Lower bound on the terminal-hitting stochastic reach
%                          avoid problem
%   optimal_input_vector - Optimal open-loop policy ((sys.input_dim) *
%                          time_horizon)-dim.  vector U = [u_0; u_1; ...; u_N]
%                          (column vector)
%
% See also computeFtLowerBoundStochReachAvoid,
% computeCcLowerBoundStochReachAvoidIterRisk,
% computeCcLowerBoundStochReachAvoidPwlRisk.
%
% Notes:
% * Delegated input handling
% * MATLAB DEPENDENCY : Uses MATLAB's Statistics and Machine Learning Toolbox
%                       Needs normpdf, normcdf, norminv for Genz's algorithm
% * EXTERNAL DEPENDENCY: Uses MPT3 and CVX
%                        Needs MPT3 for defining a controlled system and the
%                        definition of the safe and the target (polytopic) sets
% * Method 'genzps' has the following dependencies
%       * MATLAB DEPENDENCY: Uses MATLAB's Global Optimization Toolbox
%                            Needs patternsearch for gradient-free optimization
%       * EXTERNAL DEPENDENCY: Uses CVX (optional)
%                              Needs CVX to setup a convex optimization problem
%                              that initializes the patternsearch-based
%                              optimization. If CVX is unavailable, the user may
%                              provide a guess for the initialization.
%       * Specify both desired_accuracy and PSoptions or neither to use the
%         defaults 
%       * Specify an optional guess_optimal_input_vector to skip the use of CVX
% * Method 'cccpwl' has the following dependencies
%       * EXTERNAL DEPENDENCY: Uses CVX 
%                              Needs CVX to setup the convex chance-constrained
%                              problems
%       * Specify a desired_accuracy if required. Else, a default value of 1e-3
%           is used.
% * Method 'ccciter' has the following dependencies
%       * EXTERNAL DEPENDENCY: Uses CVX 
%                              Needs CVX to setup the convex chance-constrained
%                              problems
%       * Specify a desired_accuracy if required. Else, a default value of 1e-3
%           is used.
% * See @LtiSystem/getConcatMats for more information about the
%   notation used.
% * If an open_loop policy is desired arranged in increasing time columnwise,
%   use the following command
%       optimal_open_loop_control_policy = reshape(...
%           optimal_input_vector, ...
%           sys.input_dim, ...  
%           time_horizon);
% 
% =============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
%
%

    % Target tubes has polyhedra T_0, T_1, ..., T_{time_horizon}
    time_horizon = length(target_tube)-1;

    % Get half space representation of the target tube and time horizon
    % skipping the first time step
    [concat_target_tube_A, concat_target_tube_b] =...
        target_tube.concat([2 time_horizon+1]);
    n_ineq_init_set = size(target_tube(1).H,1);

    % Check if safe set contains the initial state
    if ~target_tube(1).contains(initial_state)
        % Stochastic reach-avoid probability is zero and no admissible open-loop
        % policy exists, if given an unsafe initial state
        lb_stoch_reach_avoid = 0;
        optimal_input_vector = nan(sys.input_dim * time_horizon, 1);
    else
        % Stochastic reach-avoid probability may be non-trivial
        % EXTERNAL DEPENDENCY CHECK
        if exist('cvx_begin', 'file') ~= 2
            % assert(exist('cvx_begin','file')==2, ...
            %        'SReachTools:setup_error', ...
            %        'This function uses CVX. Please get it from http://cvxr.com.');
            exc = SrtSetupError(['This function uses CVX. Please get it ', ...
                'from http://cvxr.com.']);
            throw(exc);
        end
        
        if exist('patternsearch', 'file') ~= 2
            % assert(exist('patternsearch','file')==2, ...
            %        'SReachTools:setup_error', ...
            %        'This function needs MATLAB''s Global Optimization Toolbox.');
            exc = SrtSetupError(['This function needs MATLAB''s Global ', ...
                'Optimization Toolbox.']);
            throw(exc);
        end

        if exist('normcdf', 'file') ~= 2
            % assert(exist('normcdf','file')==2, ...
            %        'SReachTools:setup_error', ...
            %        ['This function needs MATLAB''s Statistics and Machine ', ...
            %         'Learning Toolbox.']);
            exc = SrtSetupError(['This function needs MATLAB''s ', ...
                'Statistics and Machine Learning Toolbox.']);
            throw(exc);
        end

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
        [H, mean_X_sans_input, cov_X_sans_input, ~, G] = ...
            getHmatMeanCovForXSansInput(sys, ...
                initial_state, ...
                time_horizon);

        switch(lower(method))
            case 'genzps'
                % Parsing the optional arguments 
                if length(varargin) == 3
                    % First optional argument is the guess_optimal_input_vector
                    if ~isempty(varargin{1})
                        assert(isvector(varargin{1}) &&...
                               length(varargin{1}) == ...
                               sys.input_dim * time_horizon, ...
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
                               sys.input_dim * time_horizon, ...
                               'SReachTools:invalidArgs', ...
                               ['Expected a well-dimensioned row vector ', ...
                                'guess_optimal_input_vector']);
                        guess_optimal_input_vector = varargin{1};
                    else
                        guess_optimal_input_vector = [];
                    end
                    desired_accuracy = 1e-3;
                    PSoptions = psoptimset('Display', 'off');
                elseif isempty(varargin)
                    guess_optimal_input_vector = [];
                    desired_accuracy = 1e-3;
                    PSoptions = psoptimset('Display', 'off');
                else
                    exc = SrtInvalidArgsError([...
                        'guess_optimal_input_vector, desired_accuracy, ', ...
                        'and PSoptions are the only additional options.']);
                    throw(exc);
                end

                % Patternsearch and Fourier transform-based open-loop
                % underapproximation of the stochastic reach-avoid problem
                [lb_stoch_reach_avoid, optimal_input_vector] = ...
                    computeFtLowerBoundStochReachAvoid(...
                        sys, ...
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
            case 'cccpwl'
                if length(varargin) == 1
                    % Optional argument is the desired_accuracy
                    assert(isscalar(varargin{1}), ...
                           'SReachTools:invalidArgs', ...
                           'Expected a scalar value for desired_accuracy');
                    desired_accuracy = varargin{1};
                elseif isempty(varargin)
                    desired_accuracy = 1e-3;
                else
                    exc = SrtInvalidArgsError(['More inputs provided to ',...
                           'getLowerBoundStochasticReachAvoid than expected.']);
                    throw(exc);
                end
                [lb_stoch_reach_avoid, optimal_input_vector] = ...
                        computeCcLowerBoundStochReachAvoidPwlRisk(...
                            sys,...
                            time_horizon,...
                            concat_input_space_A, ... 
                            concat_input_space_b, ...
                            concat_target_tube_A, ... 
                            concat_target_tube_b, ...
                            H, ...
                            mean_X_sans_input, ...
                            cov_X_sans_input, ...
                            desired_accuracy);
            case 'ccciter'
                if length(varargin) == 1
                    % Optional argument is the desired_accuracy
                    if ~isscalar(varargin{1})
                        exc = SrtInvalidArgsError(['Expected a scalar ', ...
                            'value for desired_accuracy']);
                        throw(exc);
                    end

                    desired_accuracy = varargin{1};
                elseif isempty(varargin)
                    desired_accuracy = 1e-3;
                else
                    exc = SrtInvalidArgsError(['More inputs provided to ',...
                       'getLowerBoundStochasticReachAvoid than expected.']);
                    throw(exc);
                end
                [lb_stoch_reach_avoid, optimal_input_vector] = ...
                        computeCcLowerBoundStochReachAvoidIterRisk(...
                            sys,...
                            time_horizon,...
                            concat_input_space_A, ... 
                            concat_input_space_b, ...
                            concat_target_tube_A, ... 
                            concat_target_tube_b, ...
                            H, ...
                            mean_X_sans_input, ...
                            cov_X_sans_input, ...
                            desired_accuracy);
            case 'cccpwl-closed'
                cov_concat_disturb =kron(eye(time_horizon), ...
                                    sys.dist.parameters.covariance);
                % Direct copy-paste from cccpwl
                if length(varargin) == 1
                    % Optional argument is the desired_accuracy
                    assert(isscalar(varargin{1}), ...
                           'SReachTools:invalidArgs', ...
                           'Expected a scalar value for desired_accuracy');
                    desired_accuracy = varargin{1};
                elseif isempty(varargin)
                    desired_accuracy = 1e-3;
                else
                    exc = SrtInvalidArgsError(['More inputs provided to ',...
                           'getLowerBoundStochasticReachAvoid than expected.']);
                    throw(exc);
                end
                [lb_stoch_reach_avoid, optimal_input_vector, optimal_input_gain] = ...
                        computeClosedLoopCcLowerBoundStochReachAvoidPwlRisk(...
                            sys,...
                            time_horizon,...
                            concat_input_space_A, ... 
                            concat_input_space_b, ...
                            concat_target_tube_A, ... 
                            concat_target_tube_b, ...
                            H, ...
                            G, ...
                            mean_X_sans_input, ...
                            cov_X_sans_input, ...
                            cov_concat_disturb,...
                            desired_accuracy);
                varargout{1} = optimal_input_gain;         
            otherwise
                exc = SrtInvalidArgsError(['Unsupported method ', ...
                    'requested in getLowerBoundStochasticReachAvoid.']);
                throw(exc);
        end
        if lb_stoch_reach_avoid<0
            lb_stoch_reach_avoid = 0;
            if strcmpi(method,'ccciter') || strcmpi(method,'cccpwl')
                % Need to add sprintf (despite MATLAB's editor warning) to
                % have a new line
                % TODO: Once MILP is encoded into cccpwl, remove this.
                warning(sprintf(['Requested method returned a trivial lower',...
                    ' bound.\n Methods ''ccciter'' or ''cccpwl'' works only',...
                    ' if the maximal reach probability is above 0.5.']));
            else
                exc = SrtInternalError(['%s resulted in an infeasible ', ...
                    'optimization problem'], method);
                throw(exc);
            end
        end
    end
end
