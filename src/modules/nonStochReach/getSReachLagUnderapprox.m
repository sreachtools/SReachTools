function varargout = getSReachLagUnderapprox(sys, target_tube,...
    disturbance_set, options)
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
% [underapprox_set, underapprox_tube] = getSReachLagUnderapprox(sys,...
%       target_tube, disturbance_set)
%
% Inputs:
% -------
%   sys              - LtiSystem object
%   target_tube      - Tube object
%   disturbance_set  - Polyhedron/SReachEllipsoid object (bounded set) OR a
%                       collection of these objects which individually satisfy 
%                       the probability bound(a convex hull of the individual 
%                       results taken posteriori)
%   options          - Struct of reach set options, see SReachSetOptions
%
% Outputs:
% --------
%   overapprox_set   - Polyhedron object for the underapproximation of the 
%                      stochastic reach set
%   underapprox_tube - [Optional] Tube comprising of an underapproximation of
%                      the stochastic reach sets across the time horizon
%
% Notes:
% * From computational geometry, intersections and Minkowski differences are
%   best performed in facet representation and Minkowski sums are best
%   performed in vertex representation. However, since in this computation,
%   all three operations are required, scalability of the algorithm is severly
%   hampered, despite theoretical elegance.
% * Since box and random approaches in SReachSetOptions produce Polyhedron
%   objects for disturbance sets, we rely on MPT for all the set operations.
%   This means we do have scalability issues mentioned above.
% * For ellipsoid approach in SReachSetOptions, we seek a purely facet-based
%   operation and utilize the ray-shooting algorithm to compute a facet-based
%   underapproximation of the Minkowski sum step (via vertex-based
%   underapproximation, followed by projection, followed by convex hull
%   operation)
% * Use spreadPointsOnUnitSphere.m to compute equi_dir_vecs. 
% * equi_dir_vecs is automatically generated as part of SReachSetOptions.
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
        % minkSumInner TODO                    
        % Ainv_minus_bu = inverted_state_matrix * minus_bu;
    end
    
    if tube_length > 1
        if sys.state_dim > 2 && n_disturbances > 1
            % TODO [[url once note has been added to the google group]].'])
            warning('SReachTools:runtime', ['The convex hull operation may', ...
                'produce inconsistent or inaccurate results for systems ', ...
                'with dimensions greater than 2.'])
        end
        if options.verbose
            fprintf('Time_horizon: %d\n', tube_length-1);
        end
            
        % iterate backwards
        for itt = tube_length-1:-1:1
            % Computing effective target tube for current_time
            current_time = itt - 1;
            if options.verbose
                fprintf('Computation for time step: %d\n', current_time);
            end
            if sys.isltv()
                % Overwrite the following parameters with their
                % time-varying counterparts
                inverted_state_matrix = inv(sys.state_mat(current_time));
                minus_bu = (-sys.input_mat(current_time)) * sys.input_space;
                dist_mat = sys.dist_mat(current_time);
                % minkSumInner TODO                    
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
                    % minkSumInner TODO                    
                    % one_step_backward_reach_set = minkSumInner(...
                    %     inverted_state_matrix * new_target, Ainv_minus_bu);
                    % Use ray-shooting algorithm to underapproximate one-step
                    % backward reach set
                    effective_target = safeOneStepBackReachSet(sys,...
                        new_target, target_tube(itt), options.equi_dir_vecs,...
                        options.verbose);
                    if options.verbose
                        title(sprintf(['Safe one-step back reach set at',... 
                            ' t=%d'],current_time));
                        drawnow;
                    end
                else                                   % MPT's Polyhedron object
                    if effective_dist.isEmptySet
                        % No requirement of robustness
                        new_target = effective_target;
                    else
                        % Compute a new target set for this iteration that
                        % is robust to the disturbance
                        new_target = effective_target - effective_dist;
                    end
                    % One-step backward reach set via MPT
                    one_step_backward_reach_set = inverted_state_matrix *...
                        (new_target + minus_bu);                    

                    % Guarantee staying within target_tube by intersection
                    effective_target = intersect(one_step_backward_reach_set,...
                        target_tube(itt));
                end

                % Collect the vertices of the effective_target for each 
                % disturbance set to compute the convex hull. However, don't 
                % trigger conversion unless you really have to
                if n_disturbances > 1
                    vertices = [vertices; effective_target.V];                
                end
            end

            if n_disturbances > 1
                % Compute the polyhedron from the vertices
                effective_target_tube(itt) = Polyhedron(vertices);
            else
                effective_target_tube(itt) = effective_target;
            end
        end
    end     
    varargout{1} = effective_target_tube(1);
    if nargout > 1
        varargout{2} = effective_target_tube;
    end
end

function one_step_back_reach_polytope_underapprox = safeOneStepBackReachSet( ...
    sys, target_set, current_safe_set, equi_dir_vecs, verbose)
    % Compute the one-step back reach of target_set = {z\in X: Fz<=g} by 
    % {(x,u) \in X x U: F*(Ax + Bu)<=g}, and then projecting this set to X.
    % Here, X is the state space and U is the input space

    % Get size of equi_dir_vecs
    [dir_vecs_dim, n_vertices] = size(equi_dir_vecs);

    % Compute the polytope in X*U such that there is some x and u in the current
    % time which on application of dynamics will lie in target_set. We also add
    % the constraints defined by the current safe set as well as the input
    % space. This ensures that when we project 1) we will lie in a subset of
    % current_safe_set (effect of intersection after projection), and 2) the
    % input "selected" during the projection remains feasible.
    x_u_reaches_target_set_A = [target_set.A * [sys.state_mat sys.input_mat];
        [zeros(size(sys.input_space.A,1), sys.state_dim), sys.input_space.A];
        [current_safe_set.A zeros(size(current_safe_set.A,1), sys.input_dim)]];
    x_u_reaches_target_set_Ae = target_set.Ae * [sys.state_mat sys.input_mat];
    x_u_reaches_target_set_b = [target_set.b;
                                sys.input_space.b;
                                current_safe_set.b];
    x_u_reaches_target_set_be = target_set.be;
    x_u_reaches_target_set = Polyhedron('A', x_u_reaches_target_set_A,...
        'b', x_u_reaches_target_set_b, 'Ae', x_u_reaches_target_set_Ae,...
        'be', x_u_reaches_target_set_be);
    
    %% Find a point within the 'deep enough'
    if ~isEmptySet(x_u_reaches_target_set)
        % OPTION 1: Compute the chebyshev-center of x_u_reaches_target_set
        % center_point = x_u_reaches_target_set.chebyCenter().x;

        % OPTION 2: Compute the analytic-center of x_u_reaches_target_set
        % http://web.cvxr.com/cvx/examples/cvxbook/Ch08_geometric_probs/html/analytic_center.html
        % See Boyd and Vanderberghe's convex optimization textbook, section 8.5.3
        % cvx_begin quiet
        %     cvx_solver SDPT3;     % Require SDPT3 for exponential cone solving
        %     variable analytic_center(dir_vecs_dim,1);
        %     
        %     maximize sum(log(x_u_reaches_target_set_b -...
        %         x_u_reaches_target_set_A * analytic_center))
        %     subject to
        %         x_u_reaches_target_set_Ae * analytic_center ==...
        %             x_u_reaches_target_set_be;
        % cvx_end
        % switch cvx_status
        %     case 'Solved'
        %     case 'Inaccurate/Solved'
        %         warning('SReachTools:runTime', ['CVX returned ',...
        %             'Inaccurate/Solved, while solving a subproblem for ',...
        %             'Lagrangian underapproximation (analytic center ',...
        %             'computation). Continuing nevertheless!']);
        %     otherwise
        %         throw(SrtDevError(sprintf(['Analytic center computation',...
        %             ' failed! CVX_status: %s'], cvx_status)));
        % end
        % center_point = analytic_center;

        % OPTION 3: Compute the maximum volume inscribed ellipsoid
        % Given the polytope {a_i^T x <= b_i, i = 1,...,m }, compute an
        % ellipsoid { Bu + d | || u || <= 1 } that sits within the polytope as
        % well as has the maximum volume
        cvx_begin quiet
            cvx_solver SDPT3;
            variable mve_shape(dir_vecs_dim,dir_vecs_dim) symmetric
            variable mve_center(dir_vecs_dim)
            maximize( det_rootn( mve_shape ) )
            subject to
               for i = 1:size(x_u_reaches_target_set_A,1)
                   norm(mve_shape*x_u_reaches_target_set_A(i,:)', 2 ) +...
                        x_u_reaches_target_set_A(i,:) * mve_center <=...
                            x_u_reaches_target_set_b(i);
               end
        cvx_end
        center_point = mve_center;
        equi_dir_vecs_transf = mve_shape * equi_dir_vecs;
        equi_dir_vecs = (equi_dir_vecs_transf)./norms(equi_dir_vecs_transf,2);
    else
        throw(SrtInvalidArgsError(['Recursion led to an empty target set!',...
            ' Increasing the number of underapproximative vertices might ',...
            'help.']));
    end
    
    %% Compute vrep-based underapproximation of one-step backward reach set
    % Use the ray-shooting algorithm to compute a vertex-based
    % underapproximation of x_u_reaches_target_set
    boundary_point_mat = zeros(dir_vecs_dim, n_vertices);
    for dir_indx = 1:n_vertices
        cvx_begin quiet
            variable theta;
            variable boundary_point(dir_vecs_dim,1);            
            maximize theta;
            subject to
                theta >= 0;
                boundary_point == center_point + ...
                    theta * equi_dir_vecs(:, dir_indx);
                x_u_reaches_target_set_A * boundary_point <=...
                    x_u_reaches_target_set_b;
                x_u_reaches_target_set_Ae * boundary_point ==...
                    x_u_reaches_target_set_be;
        cvx_end        
        switch cvx_status
            case 'Solved'
                boundary_point_mat(:, dir_indx) = boundary_point;
            case 'Inaccurate/Solved'
                warning('SReachTools:runTime', ['CVX returned ',...
                    'Inaccurate/Solved, while solving a subproblem for ',...
                    'Lagrangian underapproximation. Continuing nevertheless!']);
                boundary_point_mat(:, dir_indx) = boundary_point;
            otherwise
                throw(SrtDevError(sprintf(['Underapproximation failed! ',...
                    'CVX_status: %s'], cvx_status)));
        end
    end
    before_projection_vrep_underapprox = Polyhedron('V',boundary_point_mat');
    % Compute projection onto the x space
    one_step_back_reach_polytope_underapprox =...
        before_projection_vrep_underapprox.projection(1:sys.state_dim,'vrep');    
    % Compute the half-space representation using MPT3
    one_step_back_reach_polytope_underapprox.minHRep();
    if verbose
        figure();
        if x_u_reaches_target_set.Dim <= 3
            plot(x_u_reaches_target_set,'alpha',0.3);
            hold on;plot(before_projection_vrep_underapprox);
        else
            plot(before_projection_vrep_underapprox.projection(1:3));
        end
    end    
end
