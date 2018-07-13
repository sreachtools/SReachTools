function robust_eff_target = getRobustEffTarget(sys, ...
                                                target_tube, ...
                                                disturbance, ...
                                                varargin)
% SReachTools/getRobustEffTarget Get robust Effective Target Set
% =========================================================================
%
% This function will compute the augmented effect target via the algorithm in
% the paper:
%      [[Will fill out this once paper is actually submitted]]
%
% Usage: See examples/lagrangianApproximations.m
%
% =========================================================================
%
% robust_eff_target = getRobustEffTarget(sys, ...
%                                        target_tube, ...
%                                        disturbance, ...
%                                        Name, Value)
% Inputs:
% -------
%   sys          - LtiSystem object
%   target_tube  - Cell array of Polyhedron objects 
%   disturbance  - Polyhedron object (bounded disturbance set)
% 
%   Name       | Value
%   ----------------------------------------
%   style      | 'standard', 'vrep'
%
% Outputs:
% --------
%   robust_eff_target - Polyhedron object
%
% Notes:
% * From computational geometry, intersections and Minkowski differences are
%   best performed in facet representation and Minkowski sums are best
%   performed in vertex representation. However, since in this computation,
%   all three operations are required, scalability of the algorithm is severly
%   hampered, despite theoretical elgance.
%   
% =========================================================================
% 
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 
% 

    % validate the inputs
    inpar = inputParser();
    inpar.addRequired('sys', @(x) validateattributes(x, ...
        {'LtiSystem', 'LtvSystem'}, {'nonempty'}));
    inpar.addRequired('target_tube', @(x) validateattributes(x, ...
        {'TargetTube'}, {'nonempty'}));
    inpar.addRequired('disturbance', @(x) validateattributes(x, ...
        {'Polyhedron'}, {'nonempty'}));
    inpar.addParameter('style', 'standard', @(x) validatestring(x, ...
        {'standard', 'vrep'}));

    inpar.parse(sys, target_tube, disturbance, varargin{:});

    if sys.state_dim > 4
        warning(['Because both vertex and facet representation of ', ...
            'polyhedra aer required for the necessary set recursion ', ...
            'operations computing for systems greater than 4 dimensions', ...
            'can take significant time and computational effort because ', ...
            'of the need to solve the vertex-facet enumeration problem.'])
    end
    
    % validate that all elements of the target_tube are polyhedron
    for i = 1:length(target_tube)
        validateattributes(target_tube(i), {'Polyhedron'}, {'nonempty'});
    end
    
    horizon_length = length(target_tube);
    n_disturbances = length(disturbance);

    % MPT does not support A \ P or P / A (P is Polyhedron and A is matrix)
    % so we must invert prior
    if sys.islti()
        inverted_state_matrix = inv(sys.state_mat);
        minus_bu = - sys.input_mat * sys.input_space;
    end

    effective_target_temp = target_tube(horizon_length);
    if horizon_length > 1
        if n_disturbances > 1
            if sys.state_dim > 2
                % warning(['The convex hull operation may produce ', ...
                %     'inconsistent or inaccurate results for systems with ', ...
                %     'dimensions greater than 2. See [[url once note has ', ...
                %     'been added to the google group]].'])
                warning(['The convex hull operation may produce ', ...
                    'inconsistent or inaccurate results for systems with ', ...
                    'dimensions greater than 2.'])
            end
            effective_target_tube = target_tube;
            for i = horizon_length-1:-1:1
                if ~sys.islti()
                    inverted_state_matrix = inv(sys.state_mat(i))
                    minus_bu = sys.input_mat(i) * sys.input_space;
                end

                eff_target_stopwatch = tic;

                vertices = [];
                for j = 1: n_disturbances
                    effective_dist = sys.dist_mat(i) * disturbance{j};
                    single_dist_stopwatch = tic;
                    effective_target = computeRobusteEffTargetRecursion(...
                        effective_target_tube(i+1), ...
                        effective_target_tube(i), ...
                        minus_bu, ...
                        inverted_state_matrix, ....
                        effective_dist, ...
                        inpar.Results.style);
                    
                    single_dist_comptime = toc(single_dist_stopwatch);

                    vertices = [vertices; effective_target.V];
                end
                
                effective_target_tube(i) = Polyhedron(vertices);
                
                eff_target_comptime = toc(eff_target_stopwatch);                
            end
            effective_target_temp = effective_target_tube(1);
        else
            if sys.islti()
                effective_dist = sys.dist_mat * disturbance;
            end

            for i = horizon_length-1:-1:1
                if ~sys.islti()
                    inverted_state_matrix = inv(sys.state_mat(i));
                    minus_bu = sys.input_mat(i) * sys.input_space;
                    effective_dist = sys.dist_mat(i) * disturbance;
                end
                
                eff_target_stopwatch = tic;
                effective_target_temp = ...
                    computeRobusteEffTargetRecursion(...
                        effective_target_temp, ...
                        target_tube(i), ...
                        minus_bu, ...
                        inverted_state_matrix, ...
                        effective_dist, ...
                        inpar.Results.style);
                
                eff_target_comptime = toc(eff_target_stopwatch);
            end 
        end
    end
    robust_eff_target = effective_target_temp;
end

function back_recursion_set = computeRobusteEffTargetRecursion(...
    effective_target, ...
    target_tube_set, ...
    minus_bu, ...
    inverted_state_matrix, ...
    effective_dist, ...
    style)
% SReachTools/computeRobusteEffTargetRecursion  Do the 
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
% back_recursion_set = computeRobusteEffTargetRecursion(...
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
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
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
