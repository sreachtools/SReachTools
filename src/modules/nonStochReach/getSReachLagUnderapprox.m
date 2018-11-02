function [underapprox_set, varargout] = getSReachLagUnderapprox(sys, ...
    target_tube, scaled_disturbance)
% Get underapproximation of stochastic reach set
% =========================================================================
%
% This function will compute the underapproximation of the stochastic reach
% set via Algorithm 1 in
% 
%      J. D. Gleason, A. P. Vinod, and M. M. K. Oishi. 2018. Lagrangian 
%      Approximations for Stochastic Reachability of a Target Tube. 
%      online. (2018). https://arxiv.org/abs/1810.07118
%
% Usage: See examples/lagrangianApproximations.m
%
% =========================================================================
%
% underapprox_set = getRobustEffTarget(sys, ...
%                                      target_tube, ...
%                                      disturbance, ...
%                                      Name, Value)
% Inputs:
% -------
%   sys          - LtiSystem object
%   target_tube  - Tube object
%   disturbance  - Polyhedron object (bounded disturbance set)
%
% Outputs:
% --------
%   underapprox_set - Polyhedron object
%
% Notes:
% * From computational geometry, intersections and Minkowski differences are
%   best performed in facet representation and Minkowski sums are best
%   performed in vertex representation. However, since in this computation,
%   all three operations are required, scalability of the algorithm is severly
%   hampered, despite theoretical elegance.
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
        {'Tube'}, {'nonempty'}));
    inpar.addRequired('disturbance', @(x) validateattributes(x, ...
        {'Polyhedron','SReachEllipsoid'}, {'nonempty'}));
    
    try
        inpar.parse(sys, target_tube, scaled_disturbance);
    catch cause_exc
        exc = SrtInvalidArgsError.withFunctionName();
        exc = addCause(exc, cause_exc);
        throwAsCaller(exc);
    end
    
    if sys.state_dim > 4
        warning('SReachTools:runtime',['Because both vertex and facet ', ...
            'representation of polyhedra are required for the necessary set', ...
            ' recursion operations, computing for systems greater than 4 ', ...
            'dimensions can take significant computational time and effort.']);
    end
    
    tube_length = length(target_tube);
    n_disturbances = length(scaled_disturbance);


    % initialize polyhedron array
    effective_target_tube(tube_length) = target_tube(end);
    effective_target = effective_target_tube(end);

    if sys.islti()
        % LTI System recursion
        % --------------------------------------------------------------------
        % For LTI systems we can improve the speed some by precomputing the
        % inverse of A and -B*U because they will remain consistent across
        % time

        % MPT does not support A \ P or P / A (P is Polyhedron and A is matrix)
        % so we invert prior to improve speed
        inverted_state_matrix = inv(sys.state_mat);
        minus_bu = - sys.input_mat * sys.input_space;

        if tube_length > 1
            if sys.state_dim > 2 && n_disturbances > 1
                % TODO [[url once note has been added to the google group]].'])
                warning('SReachTools:runtime', ['The convex hull operation may', ...
                    'produce inconsistent or inaccurate results for systems ', ...
                    'with dimensions greater than 2.'])
            end

            % iterate backwards
            for itt = tube_length-1:-1:1
                % Computing effective target tube for current_time
                current_time = itt - 1;

                vertices = [];
                for idist = 1:n_disturbances
                    if n_disturbances > 1
                        % Disturbance matrix has already been accounted for
                        effective_dist = scaled_disturbance{idist};
                    else
                        effective_dist = scaled_disturbance;
                    end

                    if isa(effective_dist, 'SReachEllipsoid')
                        %effective_target.minHRep();
                        new_target_A = effective_target.A;            
                        % support function of the effective_dist - vectorized to handle
                        % Implementation of Kolmanovsky's 1998-based minkowski difference
                        new_target_b = effective_target.b - ...
                            effective_dist.support_fun(new_target_A);
                        new_target = Polyhedron('H',[new_target_A new_target_b]);
                    else %if isa(effective_dist, 'Polyhedron') 
                        if effective_dist.isEmptySet
                            % No requirement of robustness
                            new_target = effective_target;
                        else
                            % Compute a new target set for this iteration that is robust to 
                            % the disturbance
                            new_target = effective_target - effective_dist;
                        end
                    end

                    % One-step backward reach set
%                     one_step_backward_reach_set = inverted_state_matrix * ...
%                         (new_target + minus_bu);
                    
                    plusInner = minkSumInner(new_target, minus_bu);
                    one_step_backward_reach_set=inverted_state_matrix*plusInner;

                    % Guarantee staying within target_tube by intersection
                    effective_target = intersect(one_step_backward_reach_set,...
                        target_tube(itt));

                    if n_disturbances > 1
                        % Don't trigger conversion unless you really have to
                        vertices = [vertices; effective_target.V];                
                    end
                end

                if n_disturbances > 1
                    effective_target_tube(itt) = Polyhedron(vertices);
                else
                    effective_target_tube(itt) = effective_target;
                end
                
%                 figure(itt);
%                 clf
%                 plot(target_tube(1).slice([3,4],zeros(2,1)),'color','y');
%                 hold on;
%                 plot(target_tube(end).slice([3,4],zeros(2,1)),'color','k');
%                 clr_string = {'w','b','r','k','g'};
%                 plot(one_step_backward_reach_set.projection([1,2]),'color',clr_string{itt},'alpha',0.2);
% %                 plot(effective_target.projection([1,2]),'color',clr_string{itt},'alpha',0.2);
%                 drawnow;              

            end
        end     
    else
        % LTV System Recursion
        % -------------------------------------------------------------------
        % The recursion is very similar to the LTI system recursion but many
        % of the sets cannot be precomputed for speed because the system
        % matrices are time-dependent

        throw(SrtInternalError('LtvSystem development is on going!'));
        % Need to fix the SReachSetLagBset generation for LtvSystem



        if tube_length > 1
            if sys.state_dim > 2 && n_disturbances > 1
                % TODO [[url once note has been added to the google group]].'])
                warning('SReachTools:runtime', ['The convex hull operation may', ...
                    'produce inconsistent or inaccurate results for systems ', ...
                    'with dimensions greater than 2.'])
            end

            % iterate backwards
            for itt = tube_length-1:-1:1
                % Computing effective target tube for current_time
                % adjustment is needed because of MATLAB's indexing at 1 instaed
                % of 0
                current_time = itt - 1;

                % run recursion for an LTI system
                % MPT does not support A \ P or P / A (P is Polyhedron and A is matrix)
                % so we must invert prior
                inverted_state_matrix = inv(sys.state_mat(current_time));
                minus_bu = - sys.input_mat(current_time) * sys.input_space;

                if iscell(scaled_disturbance)
                    for lv = 1:length(n_disturbances)
                        if isa(scaled_disturbance, 'Polyhedron');
                            scaled_disturbance{lv} = sys.dist_mat * scaled_disturbance;
                        elseif isa(scaled_disturbance, 'SReachEllipsoid')
                            scaled_disturbance{lv} = SReachEllipsoid( ...
                                sys.dist_mat * scaled_disturbance.center, ...
                                sys.dist_mat' * scaled_disturbance.shape_matrix * sys.dist_mat);
                        end
                    end
                else
                    if isa(scaled_disturbance, 'Polyhedron');
                        scaled_disturbance = sys.dist_mat * scaled_disturbance;
                    elseif isa(scaled_disturbance, 'SReachEllipsoid')
                        scaled_disturbance = SReachEllipsoid( ...
                            sys.dist_mat * scaled_disturbance.center, ...
                            sys.dist_mat' * scaled_disturbance.shape_matrix * sys.dist_mat);
                    end
                end

                vertices = [];
                for idist = 1:n_disturbances
                    if n_disturbances > 1
                        % Disturbance matrix has already been accounted for
                        effective_dist = scaled_disturbance{idist};
                    else
                        effective_dist = scaled_disturbance;
                    end

                    if isa(effective_dist, 'Polyhedron') && ....
                       effective_dist.isEmptySet
                        % No requirement of robustness
                        new_target = effective_target;
                    elseif isa(effective_dist, 'SReachEllipsoid')
                        %effective_target.minHRep();
                        new_target_A = effective_target.A;            
                        % support function of the effective_dist - vectorized to handle
                        % Implementation of Kolmanovsky's 1998-based minkowski difference
                        new_target_b = effective_target.b -...
                            effective_dist.support_fun(new_target_A);
                        new_target = Polyhedron('H',[new_target_A new_target_b]);
                    else
                        % Compute a new target set for this iteration that is robust to 
                        % the disturbance
                        new_target = effective_target - effective_dist;
                    end

                    % One-step backward reach set
                    one_step_backward_reach_set = inverted_state_matrix * ...
                        (new_target + minus_bu);


                    % Guarantee staying within target_tube by intersection
                    effective_target = intersect(...
                        one_step_backward_reach_set, ...
                        target_tube(itt));

                    if n_disturbances > 1
                        % Don't trigger conversion unless you really have to
                        vertices = [vertices; effective_target.V];                
                    end
                end

                if n_disturbances > 1
                    effective_target_tube(itt) = Polyhedron(vertices);
                else
                    effective_target_tube(itt) = effective_target;
                end
            end
        end    
    end

    underapprox_set = effective_target_tube(1);
    if tube_length > 1
        varargout{1} = effective_target_tube(2:end);
    end
end

function back_recursion_set = computeUnderapproxsetRecursion(...
    effective_target, ...
    target_tube_set, ...
    minus_bu, ...
    inverted_state_matrix, ...
    effective_dist, ...
    style)
% Do the one set backward recursion to obtain the effective target
% =========================================================================
%
% THIS FUNCTION IS NOT CURRENTLY IN USE. IT IS NOT YET DELETED BUT IS NOT CALLED
% Nested function to perform the one-step recursion of Algorithm 1 from:
% 
%      J. D. Gleason, A. P. Vinod, and M. M. K. Oishi. 2018. Lagrangian 
%      Approximations for Stochastic Reachability of a Target Tube. 
%      online. (2018). https://arxiv.org/abs/1810.07118
% 
% Usage: Nested function in getRobustEffTarget
%
% =========================================================================
% 
% back_recursion_set = computeUnderapproxsetRecursion(...
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
    
    if isa(effective_dist, 'Polyhedron') && effective_dist.isEmptySet
        % No requirement of robustness
        new_target = effective_target;
    elseif isa(effective_dist, 'SReachEllipsoid')
        %effective_target.minHRep();
        new_target_A = effective_target.A;            
        % support function of the effective_dist - vectorized to handle
        % Implementation of Kolmanovsky's 1998-based minkowski difference
        new_target_b = effective_target.b -...
            effective_dist.support_fun(new_target_A);
        new_target = Polyhedron('H',[new_target_A new_target_b]);
    else
        % Compute a new target set for this iteration that is robust to 
        % the disturbance
        new_target = effective_target - effective_dist;
    end
            
    % choose style
    switch(style)
        case 'standard'
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
            throw(SrtInternalError('Unhandled option'));
    end

end
