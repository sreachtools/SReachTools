function grid_probability = getDynProgSolForTargetTube(sys, ...
    state_grid, input_grid, target_tube, varargin)
% SReach/stochasticReachAvoid/getDynProgSolForTargetTube
% ============================================================================
%
% The function computes the probability of staying in a target tube defined
% on a particular state stace grid. The dynamic programming recursion can be
% found in 
%   
% S. Summers and J. Lygeros, "Verification of discrete time stochastic hybrid 
% systems: A stochastic reach-avoid decision problem," Automatica, vol. 46,
% no. 12, pp. 1951--1961, 2010.
%
% The problem of examining the reachability of a target tube can be found in
% a work that we intend to publish soon :D
%
% Usage: See example doubleIntegratorDynmaicProgramming.m
%
% ============================================================================
%
% grid_probability = getDynProgSolForTargetTube(sys, ...
%     state_grid, input_grid, target_tube, varargin)
% 
% Inputs:
% -------
%   sys         - LtiSystem object
%   state_grid  - SpaceGrid object
%   input_grid  - InputGrid object
%   target_tube - Cell array of Polyhedron objects
%
% Outputs:
% --------
%   grid_probability - Nx1 Array of probability values, where N is equivalent
%                      to size(state_grid, 1)
%
% Notes:
%   - WARNING: Dynamic programming suffers from the curse of dimensionality! As
%     such, this code will effective and with reasonable speed compute dynamic
%     programming solutions for 2-dimensional systems with Gaussian 
%     disturbances. However, for 3-dimensional systems or larger the required
%     computation time, and memory, will exponentially grow to the point that 
%     the simulation will take longer than it took for me to get my PhD.
%   
%   - Currently this back propagation, and subsequently the entire dynamic 
%     programming recursion, only works for Gaussian disturbances.
% 
% ============================================================================
% 
%   This function is part of the Stochastic Optimal Control Toolbox.
%   License for the use of this function is given in
%        https://github.com/abyvinod/SReach/blob/master/LICENSE
%


    if length(varargin) > 0
        options = processDynamicProgrammingOptions(varargin{:});
    else
        options = processDynamicProgrammingOptions();
    end

    % check inputs
    validateattributes(sys, {'LtiSystem'}, {'nonempty'})
    validateattributes(state_grid, {'SpaceGrid'}, {'nonempty'})
    validateattributes(input_grid, {'InputGrid'}, {'nonempty'})
    validateattributes(target_tube, {'cell'}, {'nonempty', 'vector'});
    for i = 1:length(target_tube)
        validateattributes(target_tube{i}, {'Polyhedron'}, {'nonempty'});
    end
    
    % need a stochastic disturbance of type Gaussian for probability computation
    % right now
    validateattributes(sys.disturbance, {'StochasticDisturbance'}, ...
        {'nonempty'});
    if ~strcmpi(sys.disturbance.type, 'Gaussian')
        error('SReach:invalidArgs', ['Dynamic programming methods can ', ...
            'only handle stochastic disturbances of type Gaussian']);
    end
    
    % start with initialization of probability of 1s for everything in final
    % target tube
    grid_probability = zeros(size(state_grid.grid, 1), 1);
    ind_vector = state_grid.getIndicatorVectorForSet(target_tube{end});
    grid_probability(ind_vector) = 1;
    
    n_targets = length(target_tube);
    options.timer = tic;
    options.status.total_targets = n_targets;
    if options.verbose
        fprintf(['Starting dynamic programming computation for target ', ...
            'tube...\n\n'])
    end
    if n_targets > 1
        options.timer = tic;
        for i = 1:n_targets-1
            options.status.current_target = i;
            options.status.current_grid_index = 0;
            if options.verbose
                fprintf(['Target %d/%d: - simulation time %f ', ...
                    'seconds\n'], i, n_targets-1, toc(options.timer))
            end

            % compute the 1-step back propagation
            % can eventually allow for an option that will return the grid 
            % probability after each time instant, saving some computation time
            % if you are trying to generate images for multiple time points
            grid_probability = computeDynProgBackPropagation(...
                sys, ...
                state_grid, ...
                input_grid, ...
                grid_probability, ...
                target_tube{n_targets-i}, ...
                options);
        end
    end

    if strcmpi(options.status.timer.Running, 'on')
        stop(options.status.timer);
    end

end

function options = processDynamicProgrammingOptions(varargin)
% SReach/getDynProgSolForTargetTube/processDynamicProgrammingOptions
% Process options given to the Dynamic programming solver
% ============================================================================
%
% Process options given to the getDynProgSolForTargetTube
% solver function
% 
% Currently under construction so is not used
%   
% ============================================================================
%
% options = processDynamicProgrammingOptions(varargin)
% 
% Inputs:
% -------
%   varargin - Options, either in singluar option or Name Value pairs
%
% Outputs:
% --------
%   options - Options struct
%
% Notes:
%   * Under construction, check back soon...
% 
% ============================================================================
% 
%   This function is part of the Stochastic Optimal Control Toolbox.
%   License for the use of this function is given in
%        https://github.com/abyvinod/SReach/blob/master/LICENSE
%    

    options = struct();

    % Default options
    options.timer   = [];
    options.verbose = false;
    options.status.show = false;
    options.status.period = NaN;
    options.status.current_target = 0;
    options.status.total_targets = 0;
    options.status.timer = timer();

    if ~isempty(varargin)
        i = 1;
        while i <= length(varargin)
            if strcmpi(varargin{i}, 'verbose')
                options.verbose = true;
                i = i + 1;
            elseif strcmpi(varargin{i}, 'displaystatus')
                options.status.show = true;
                options.status.period = varargin{i+1};
                options.status.timer.Period = varargin{i+1};
                options.status.timer.ExecutionMode = 'fixedRate';
                i = i + 2;
            end
        end
    end

end