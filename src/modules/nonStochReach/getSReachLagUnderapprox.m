function [underapprox_set, varargout] = getSReachLagUnderapprox(sys, ...
    target_tube, disturbance_set)
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
% underapprox_set = getSReachLagUnderapprox(sys, target_tube, disturbance_set)
%
% Inputs:
% -------
%   sys              - LtiSystem object
%   target_tube      - Tube object
%   disturbance_set  - Polyhedron/SReachEllipsoid object (bounded set) OR a
%                       collection of these objects which individually satisfy 
%                       the probability bound(a convex hull of the individual 
%                       results taken posteriori)
%
% Outputs:
% --------
%   underapprox_set  - Polyhedron object
%   effective_target_tube
%                    - [Optional] Tube comprising of an underapproximation
%                           of the stochastic reach sets across the time horizon
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
        inpar.parse(sys, target_tube, disturbance_set);
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
    n_disturbances = length(disturbance_set);


    % initialize polyhedron array
    effective_target_tube(tube_length) = target_tube(end);
    effective_target = effective_target_tube(end);

    if sys.islti()
        inverted_state_matrix = inv(sys.state_mat);
        minus_bu = (-sys.input_mat) * sys.input_space;
        dist_mat = sys.dist_mat;
        % % For Minkowski inner approximation: TODO
        % Ainv_minus_bu = inverted_state_matrix * minus_bu;
    end
    
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
            if sys.isltv()
                % Overwrite the following parameters with their
                % time-varying counterparts
                inverted_state_matrix = inv(sys.state_mat(current_time));
                minus_bu = (-sys.input_mat(current_time)) * sys.input_space;
                dist_mat = sys.dist_mat(current_time);
                % % For Minkowski inner approximation: TODO
                % Ainv_minus_bu = inverted_state_matrix * minus_bu;
            end
                
            vertices = [];
            
            for idist = 1:n_disturbances
                % Account for disturbance matrix
                if n_disturbances > 1
                    effective_dist = dist_mat * disturbance_set{idist};
                else
                    effective_dist = dist_mat * disturbance_set;
                end
                
                if isa(effective_dist, 'SReachEllipsoid')
                    % support function of the effective_dist - vectorized to
                    % handle Implementation of Kolmanovsky's 1998-based
                    % minkowski difference
                    new_target_A = effective_target.A;            
                    new_target_b = effective_target.b - ...
                        effective_dist.support_fun(new_target_A);
                    new_target= Polyhedron('H',[new_target_A new_target_b]);
                else                                   % MPT's Polyhedron object
                    if effective_dist.isEmptySet
                        % No requirement of robustness
                        new_target = effective_target;
                    else
                        % Compute a new target set for this iteration that
                        % is robust to the disturbance
                        new_target = effective_target - effective_dist;
                    end
                end

                % One-step backward reach set via MPT
                one_step_backward_reach_set = inverted_state_matrix *...
                        (new_target + minus_bu);                    
%               % ALTERNATIVELY, use minkSumInner TODO                    
%                 one_step_backward_reach_set = minkSumInner(...
%                     inverted_state_matrix * new_target, Ainv_minus_bu);

                % Guarantee staying within target_tube by intersection
                effective_target = intersect(one_step_backward_reach_set,...
                    target_tube(itt));

                % Collect the vertices of the effective_target for each 
                % disturbance set to compute the convex hull. However, don't 
                % trigger conversion unless you really have to
                if n_disturbances > 1
                    vertices = [vertices; effective_target.V];                
                end
            end

            if n_disturbances > 1
                % Compute the convex hull
                effective_target_tube(itt) = Polyhedron(vertices);
            else
                effective_target_tube(itt) = effective_target;
            end
        end
    end     

    underapprox_set = effective_target_tube(1);
    if tube_length > 1 && nargout > 1
        varargout{1} = effective_target_tube;
    end
end