function robust_eff_target = getRobustEffTarget(sys, ...
                                                target_tube, ...
                                                disturbance, ...
                                                options)
% SReach/getRobustEffTarget Get robust Effective Target Set
% =========================================================================
%
% This function will compute the augmented effect target via the algorithm in
% the paper:
%      [[Will fill out this once paper is actually submitted]]
%
% Usage:
% ------
% Will create example to demonstrate use
%
% =========================================================================
%
% Inputs:
% -------
%   sys          - LtiSystem object
%   target_tube  - Cell array of Polyhedron objects 
%   disturbance  - Polyhedron object (bounded disturbance set)
%   options      - (optional) Logging options, currently in development
%
% Outputs:
% --------
%   robust_eff_target - Polyhedron object
%
% Notes:
%   - From computational geometry, intersections and Minkowski differences are
%     best performed in facet representation and Minkowski sums are best
%     performed in vertex representation. However, since in this computation,
%     all three operations are required, scalability of the algorithm is severly
%     hampered, despite theoretical elgance.
%   
% =========================================================================
% 
%   This function is part of the Stochastic Optimal Control Toolbox.
%   License for the use of this function is given in
%        https://github.com/abyvinod/SReach/blob/master/LICENSE
% 
% 

    % validate the inputs
    if nargin < 4
        options = getDefaultOptions();
    elseif nargin < 3
        error('SReach:Internal', 'No enough input arguments.');
    end
    
    validateattributes(sys, {'LtiSystem'}, {'nonempty'});
    validateattributes(target_tube, {'cell'}, {'nonempty'});
    validateattributes(disturbance, {'Polyhedron', 'cell'}, {'nonempty'});

    if sys.state_dimension > 4
        warning(['Because both vertex and facet representation of ', ...
            'polyhedra aer required for the necessary set recursion ', ...
            'operations computing for systems greater than 4 dimensions', ...
            'can take significant time and computational effort because ', ...
            'of the need to solve the vertex-facet enumeration problem.'])
    end
    
    % validate that all elements of the target_tube are polyhedron
    for i = 1:length(target_tube)
        validateattributes(target_tube{i}, {'Polyhedron'}, {'nonempty'});
    end
    
    horizon_length = length(target_tube);
    n_disturbances = length(disturbance);

    % MPT does not support A \ P or P / A (P is Polyhedron and A is matrix)
    % so we must invert prior
    inverted_state_matrix = inv(sys.state_matrix);
    minus_bu = - sys.input_matrix * sys.input_space;
    effective_target_temp = target_tube{horizon_length};
    if horizon_length > 1
        if n_disturbances > 1
            if sys.state_dimension > 2
                warning(['The convex hull operation may produce ', ...
                    'inconsistent or inaccurate results for systems with ', ...
                    'dimensions greater than 2. See [[url once note has ', ...
                    'been added to the google group]].'])
            end
            effective_target_tube = target_tube;
            for i = horizon_length-1:-1:1
                eff_target_stopwatch = tic;
                if options.verbose == 1
                   fprintf('Computing effective target for k=%d...\n', i-1)
                end

                vertices = [];
                for j = 1: n_disturbances
                    effective_dist = sys.disturbance_matrix * disturbance{j};
                    single_dist_stopwatch = tic;
                    if options.verbose == 1
                       fprintf(['    Computing effective target for ', ...
                           'disturbance %d/%d... '], j, n_disturbances);
                    end
                    effective_target = performRobustEffectiveTargetRecursion(...
                        effective_target_tube{i+1}, ...
                        effective_target_tube{i}, ...
                        minus_bu, ...
                        inverted_state_matrix, ....
                        effective_dist, ...
                        options.style);
                    
                    single_dist_comptime = toc(single_dist_stopwatch);
                    if options.verbose == 1
                        fprintf('Computation time: %.3f\n', ...
                            single_dist_comptime);
                    end

                    vertices = [vertices; effective_target.V];
                end
                
                effective_target_tube{i} = Polyhedron(vertices);
                
                eff_target_comptime = toc(eff_target_stopwatch);
                if options.verbose == 1
                    fprintf(['Total computation time for diturbances: ', ...
                        '%.3f\n'], eff_target_comptime);
                end
                
            end
            effective_target_temp = effective_target_tube{1};
        else
            effective_dist = sys.disturbance_matrix * disturbance;
            for i = horizon_length-1:-1:1
                eff_target_stopwatch = tic;
                if options.verbose == 1
                   fprintf('Computing effective target for k=%d... ', i-1)
                end
                effective_target_temp = ...
                    performRobustEffectiveTargetRecursion(...
                        effective_target_temp, ...
                        target_tube{i}, ...
                        minus_bu, ...
                        inverted_state_matrix, ...
                        effective_dist, ...
                        options.style);
                
                eff_target_comptime = toc(eff_target_stopwatch);
                if options.verbose == 1
                    fprintf('Computation time: %.3f\n', eff_target_comptime);
                end
            end 
        end
    end
    robust_eff_target = effective_target_temp;
end

function back_recursion_set = performRobustEffectiveTargetRecursion(...
    effective_target, ...
    target_tube_set, ...
    minus_bu, ...
    inverted_state_matrix, ...
    effective_dist, ...
    style)
% SReach/performRobustEffectiveTargetRecursion  Do the 
% one set backward recursion to obtain the effective target
% =========================================================================
%
% Nested function to perform the one-step recursion for obtaining the robust
% effective target set. See recursion from
%   [[Will fill out this once paper is actually submitted]]
% 
% Usage
% -----
% Nested function in getRobustEffTarget
%
% =========================================================================
% 
% back_recursion_set = performRobustEffectiveTargetRecursion(...
%     effective_target, ...
%     target_tube_set, ...
%     minus_bu, ...
%     inverted_state_matrix, ...
%     effective_dist, ...
%     style);
% 
% Inputs
% ------
%   effective_target      - Current effective target (polyhedron)    
%   target_tube_set       - k-time target tube set (polyhedron)
%   minus_bu              - Set (-input matrix * input space) (polyhedron)
%   inverted_state_matrix - Inversion of state matrix A^{-1}
%   effective_dist        - Effective disturbance (disturbance matrix * 
%                           bounded disturbance set) (polyhedron)
%   style                 - Choice of recurcion excecution style (string)
% 
%   Acceptable styles:
%       * 'standard' - Standard MPT execution
%       * 'vrep'     - Computes most operations in vertex representation
% 
% Outputs
% -------
% back_recursion_set - Backward recursion result (polyhedron)
% 
% =========================================================================
% 
%   This function is part of the Stochastic Optimal Control Toolbox.
%   License for the use of this function is given in
%        https://github.com/abyvinod/SReach/blob/master/LICENSE
% 
% 
    
    % choose style
    switch(style)
        case 'standard'
            if effective_dist.isEmptySet
                % No requirement of robustness
                new_target = effective_target;
            else
                % Compute a new target set for this iteration that is robust to 
                % the disturbance
                new_target = effective_target - ...
                    effective_dist;
            end

            % One-step backward reach set
            one_step_backward_reach_set = inverted_state_matrix * ...
                (new_target + minus_bu);


            % Guarantee staying within target_tube by intersection
            back_recursion_set = intersect(...
                one_step_backward_reach_set, ...
                target_tube_set);

        case 'vrep'
            % running the computation using a vrep style is an augmentation that
            % tries to keep the polytope in vertex representation for as many
            % possible operations
            if effective_dist.isEmptySet
                % No requirement of robustness
                new_target = effective_target;
            else
                % Compute a new target set for this iteration that is robust to 
                % the disturbance
                new_target = effective_target - effective_dist;
            end

            % Miniminizing H before MPT makes us jump to V
            minHRep(new_target);

            % Will return a V-Polytope 
            % --- CDDMEX might cause failure 
            % here because of the vertex-facet enumeration
            temp_reach_set = plus(new_target, minus_bu, 'vrep');     

            % Will return a V-Polytope
            reach_set = affineMap(temp_reach_set, inverted_state_matrix, ...
                'vrep');

            % Miniminizing V before MPT makes us jump to H                  
            minVRep(reach_set);

            % Will return a H-Polytope 
            % --- CDDMEX might cause failure
            % here because of the vertex-facet enumeration
            back_recursion_set = intersect(reach_set, target_tube_set);
        otherwise
            assert(false, 'Unhandled option')
    end

end

function options = getDefaultOptions()
% SReach/getRobustEffTarget/getDefaultOptions  Get default solver options
% =========================================================================
%
% Nested function to obtain default solver/logging options
% 
%   Usage
%   -----
%   Nested function
%
% =========================================================================
% 
% options = GETDEFAULTOPTIONS()
% 
% Inputs
% ------
% None
% 
% Outputs
% -------
% options - options struct
% 
% =========================================================================
% 
%   This function is part of the Stochastic Optimal Control Toolbox.
%   License for the use of this function is given in
%        https://github.com/abyvinod/SReach/blob/master/LICENSE
% 
% 

options = struct();
options.verbose = 0;
options.style = 'standard';
options.suppress_warning = false;

end