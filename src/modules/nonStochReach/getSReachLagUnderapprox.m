function varargout = getSReachLagUnderapprox(sys, target_tube,dist_set, options)
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
%   dist_set         - Polyhedron/SReachEllipsoid object (bounded set) OR a
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
%        https://sreachtools.github.io/license/
% 
% 

    % validate the inputs
    inpar = inputParser();
    inpar.addRequired('sys', @(x) validateattributes(x, ...
        {'LtiSystem', 'LtvSystem'}, {'nonempty'}));
    inpar.addRequired('target_tube', @(x) validateattributes(x, ...
        {'Tube'}, {'nonempty'}));
    inpar.addRequired('dist_set', @(x) validateattributes(x, ...
        {'Polyhedron','SReachEllipsoid'}, {'nonempty'}));
    
    try
        inpar.parse(sys, target_tube, dist_set);
    catch cause_exc
        exc = SrtInvalidArgsError.withFunctionName();
        exc = addCause(exc, cause_exc);
        throwAsCaller(exc);
    end
 
    % Check if prob_str and method_str of options are consistent        
    if ~strcmpi(options.prob_str, 'term')
        throwAsCaller(...
            SrtInvalidArgsError('Mismatch in prob_str in the options'));
    end
    if ~strcmpi(options.method_str, 'lag-under')
        throwAsCaller(...
            SrtInvalidArgsError('Mismatch in method_str in the options'));
    end            
    
    if sys.state_dim > 4
        warning('SReachTools:runtime',['Because both vertex and facet ', ...
            'representation of polyhedra are required for the necessary set', ...
            ' recursion operations, computing for systems greater than 4 ', ...
            'dimensions can take significant computational time and effort.']);
    end
    
    tube_length = length(target_tube);
    n_disturbances = length(dist_set);


    % initialize polyhedron array
    effective_target_tube(tube_length) = target_tube(end);
    effective_target = effective_target_tube(end);

    if sys.islti()
        inverted_state_matrix = inv(sys.state_mat);
        minus_bu = (-sys.input_mat) * sys.input_space;
        dist_mat = sys.dist_mat;
    end
    
    if tube_length > 1
        if options.verbose >= 1
            fprintf('Time_horizon: %d\n', tube_length-1);
        end
            
        % iterate backwards
        for itt = tube_length-1:-1:1
            % Computing effective target tube for current_time
            current_time = itt - 1;
            if options.verbose >= 1
                if options.verbose >= 2
                    fprintf('\n');
                end
                fprintf('Computation for time step: %d\n', current_time);
            end
            if sys.isltv()
                % Overwrite the following parameters with their
                % time-varying counterparts
                inverted_state_matrix = inv(sys.state_mat(current_time));
                minus_bu = (-sys.input_mat(current_time)) * sys.input_space;
                dist_mat = sys.dist_mat(current_time);
            end
                
            vertices = [];
            
            for idist = 1:n_disturbances
                % Account for disturbance matrix
                if n_disturbances > 1
                    effective_dist = dist_mat * dist_set{idist};
                else
                    effective_dist = dist_mat * dist_set;
                end
                
                % support function of the effective_dist - vectorized to
                % handle Implementation of Kolmanovsky's 1998-based
                % minkowski difference
                new_target_A = effective_target.A;            
                if isa(effective_dist, 'SReachEllipsoid')
                    new_target_b = effective_target.b - ...
                        effective_dist.support(new_target_A');
                elseif isa(effective_dist, 'Polyhedron')
                    if effective_dist.isEmptySet()
                        new_target_b = effective_target.b;
                    else
                        new_target_b = effective_target.b - ...
                            effective_dist.support(new_target_A');
                    end
                else
                    throwAsCaller(SrtInvalidArgsError(sprintf(['Invalid ',...
                        'disturbance object (%s). Expected SReachEllipsoid/',...
                        'Polyhedron.'], class(effective_dist))));
                end
                new_target= Polyhedron('H',[new_target_A new_target_b]);

                switch lower(options.compute_style)
                    case 'support'
                        % Use ray-shooting algorithm to underapproximate 
                        % one-step backward reach set
                        if options.verbose >= 2
                            timerVal = tic;
                        end
                        effective_target = safeOneStepBackReachSet(itt-1, ...
                            sys, new_target, target_tube(itt), options);
                        if options.verbose >= 2
                            fprintf('Time for oneStepBack: %1.3f s\n',...
                                toc(timerVal));
                        end
                        if options.verbose == 3
                            title(sprintf(['Safe one-step back reach set at',... 
                                ' t=%d'],current_time));
                            drawnow;
                        end
                    case 'vfmethod'
                        switch options.vf_enum_method
                            case 'cdd'
                                %% One-step backward reach set via MPT
                                one_step_backward_reach_set = ...
                                    inverted_state_matrix * ...
                                        (new_target + minus_bu);                    
                                % Guarantee staying within target_tube by
                                % intersection
                                effective_target = intersect(...
                                    one_step_backward_reach_set, ...
                                    target_tube(itt));
                            case 'lrs'
                                %% One-step backward reach set via LRS
                                effective_target = lrsOneStepBackReachSet(...
                                    itt - 1, sys, new_target, target_tube(itt), ...
                                    options.verbose);
                        end
                    otherwise
                        throwAsCaller(SrtInvalidArgsError(['Invalid ',...
                            'computation style specified']));
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
    current_time, sys, target_set, current_safe_set, options)
    % Compute the one-step back reach of target_set = {z\in X: Fz<=g} by 
    % {(x,u) \in X x U: F*(Ax + Bu)<=g}, and then projecting this set to X.
    % Here, X is the state space and U is the input space

    if isempty(options.equi_dir_vecs)
        throwAsCaller(SrtInvalidArgsError(['Expected non-empty ',...
            'equi_dir_vecs. Faulty options structure provided!']));
    end
    
    % Get size of equi_dir_vecs
    [dir_vecs_dim, n_vertices] = size(options.equi_dir_vecs);
    
    if dir_vecs_dim ~= (sys.state_dim + sys.input_dim) || n_vertices < 3
        throwAsCaller(SrtInvalidArgsError(['Expected (sys.state_dim + ',...
            'sys.input_dim)-dimensional collection of column vectors. ',...
            'Faulty options structure provided!']));        
    end
    
    if options.verbose >= 2
        timerVal = tic;
    end
    % Compute the polytope in X*U such that there is some x and u in the current
    % time which on application of dynamics will lie in target_set. We also add
    % the constraints defined by the current safe set as well as the input
    % space. This ensures that when we project 1) we will lie in a subset of
    % current_safe_set (effect of intersection after projection), and 2) the
    % input "selected" during the projection remains feasible.
    if sys.islti
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
    else
        x_u_reaches_target_set_A = [target_set.A * ... 
                [sys.state_mat(current_time) sys.input_mat(current_time)];
            [zeros(size(sys.input_space.A,1), sys.state_dim), sys.input_space.A];
            [current_safe_set.A zeros(size(current_safe_set.A,1), sys.input_dim)]];
        x_u_reaches_target_set_Ae = target_set.Ae * ... 
                [sys.state_mat(current_time) sys.input_mat(current_time)];
        x_u_reaches_target_set_b = [target_set.b;
                                    sys.input_space.b;
                                    current_safe_set.b];
        x_u_reaches_target_set_be = target_set.be;
        x_u_reaches_target_set = Polyhedron('A', x_u_reaches_target_set_A,...
            'b', x_u_reaches_target_set_b, 'Ae', x_u_reaches_target_set_Ae,...
            'be', x_u_reaches_target_set_be);
    end
    if options.verbose >= 2
        fprintf('Time to setup x_u_reaches_target_set: %1.3f s\n', ...
            toc(timerVal));
    end
    
    %% Find a point within the 'deep enough'
    if ~isEmptySet(x_u_reaches_target_set)
        % OPTION 1: Compute the chebyshev-center of x_u_reaches_target_set
        % >>> Abandoned because the point is not a center and no transformation
        %     of the equi-spaced vectors is available
        % OPTION 2: Compute the analytic-center of x_u_reaches_target_set
        % >>> Abandoned because no transformation of the equi-spaced vectors is 
        %     available
        % OPTION 3: Compute the maximum volume inscribed ellipsoid
        % Given the polytope {a_i^T x <= b_i, i = 1,...,m }, compute an
        % ellipsoid { Bu + d | || u || <= 1 } that sits within the polytope as
        % well as has the maximum volume
        if options.verbose >= 2
            timerVal = tic;
        end
        cvx_begin quiet
            cvx_solver SDPT3;
            variable mve_shape(dir_vecs_dim,dir_vecs_dim) symmetric;
            variable mve_center(dir_vecs_dim);
            maximize( det_rootn( mve_shape ) )
            subject to
               norms(mve_shape*x_u_reaches_target_set_A', 2)' +...
                    x_u_reaches_target_set_A * mve_center <=...
                        x_u_reaches_target_set_b;
        cvx_end
        center_point = mve_center;
        equi_dir_vecs_transf = mve_shape * options.equi_dir_vecs;
        transf_equi_dir_vecs = (equi_dir_vecs_transf)./ ...
            norms(equi_dir_vecs_transf,2);
        if options.verbose >= 2
            fprintf('Time for MVE computation      : %1.3f s\n',toc(timerVal));
        end
    else
        throw(SrtRuntimeError(['Recursion led to an empty target set!',...
            ' Increasing the number of underapproximative vertices might ',...
            'help.']));
    end
    
    %% Compute vrep-based underapproximation of one-step backward reach set
    % Use the ray-shooting algorithm to compute a vertex-based
    % underapproximation of x_u_reaches_target_set
    boundary_point_mat = zeros(dir_vecs_dim, n_vertices);
    if options.verbose >= 2
        timerVal = tic;
    end
    for dir_indx = 1:n_vertices
        dir_vec = transf_equi_dir_vecs(:, dir_indx);
        
        % OPTION 1: Formulate a LP in CVX and solve it
        % >>> Abandoned because of the compute overhead
        % OPTION 2: Bisection
        boundary_point = @(theta) center_point + theta * dir_vec;
        contains_check = @(theta) all(x_u_reaches_target_set_A *...
             boundary_point(theta) <= x_u_reaches_target_set_b);
        if any(abs(x_u_reaches_target_set_Ae * (center_point + dir_vec) -...
                x_u_reaches_target_set_be) > 1e-10)
            % Skipping the direction vector since it does not satisfy the
            % equality constraint
            warning('SReachTools:runTime', ['Given vector did not satisfy ',...
                    'the equality constraint imposed. Skipping...']);
            continue;
        else
            % Bracketting the bounds for the bisection
            bisection_lb = 0;
            bisection_ub = 1;
            while contains_check(bisection_ub)
                bisection_lb = bisection_ub;
                bisection_ub = 2 * bisection_ub;                
            end
            if contains_check(bisection_lb) == 1 &&...
                    contains_check(bisection_ub) == 0
                % all ok
            else
                throw(SrtDevError('Bisection algorithm failed!'));
            end
        end

        % The threshold of 1e-4 was heuristically arrived at! CDDMEX had
        % trouble computing the facet form for vertices created with eps
        % smaller than 1e-4
        % https://github.com/cddlib/cddlib#numerical-problems
        while abs(bisection_lb-bisection_ub) > 1e-4
            bisection_test = (bisection_lb + bisection_ub) / 2;
            if contains_check(bisection_test)
                % Increase the lower bound since the given test point passed
                bisection_lb = bisection_test;
            else
                % Decrease the upper bound since the given test point failed
                bisection_ub = bisection_test;
            end
        end
        boundary_point_mat(:, dir_indx) = boundary_point(bisection_lb);
    end
    if options.verbose >= 2
        fprintf('Time to get v-rep underapprox.: %1.3f s\n',toc(timerVal));
    end
    % Compute projection onto the x space --- ignore the u dimensions and
    % define the polytope using the resulting vertices (their convex hull)
    if options.verbose >= 2
        timerVal = tic;
    end
    switch options.vf_enum_method
        case 'lrs' % LRS approach
            one_step_back_reach_polytope_underapprox_V =...
                vertexReduction(boundary_point_mat(1:sys.state_dim,:)');
        case 'cdd' % CDD approach (MPT3's native method)
            % Construct a polyhedron
            one_step_back_reach_polytope_underapprox = ...
                Polyhedron('V',boundary_point_mat(1:sys.state_dim,:)');
            % Remove the redundant vertices
            one_step_back_reach_polytope_underapprox.minVRep();   
            % Extract the irredundant vertices
            one_step_back_reach_polytope_underapprox_V =...
                one_step_back_reach_polytope_underapprox.V;
    end
    if options.verbose >= 2
        fprintf('Time to get minimal v-rep     : %1.3f s\n',toc(timerVal));
    end
    % Compute the half-space representation using MPT3
    if options.verbose >= 2
        timerVal = tic;
    end
    switch options.vf_enum_method
        case 'lrs' % LRS approach
            [one_step_back_reach_polytope_underapprox_A_red, ...
                one_step_back_reach_polytope_underapprox_b_red] = ...
                    facetEnumeration( ...
                        one_step_back_reach_polytope_underapprox_V, ...
                        ones(size( ...
                          one_step_back_reach_polytope_underapprox_V,1),1));
            [one_step_back_reach_polytope_underapprox_A, ...
                one_step_back_reach_polytope_underapprox_b] = ...
                    inequalityReduction(...
                        one_step_back_reach_polytope_underapprox_A_red ,...
                        one_step_back_reach_polytope_underapprox_b_red);
        case 'cdd' % CDD approach (MPT3's native method)
            try
                % Construct the irredundant half-space representation
                one_step_back_reach_polytope_underapprox.minHRep();
                % Extract the half-space representation
                one_step_back_reach_polytope_underapprox_A =...
                    one_step_back_reach_polytope_underapprox.A;
                one_step_back_reach_polytope_underapprox_b =...
                    one_step_back_reach_polytope_underapprox.b;
            catch ME
                if strcmpi(ME.identifier, 'SReachTools:badConv')
                    throw(SrtDevError(['Conversion from vertex to facet ',...
                    'worked, but the halfspace returned was garbage!\nThis',...
                    ' is most likely due to numerical issues in ',...
                    'convhulln/CDDMEX in MPT.']));
                else
                    throw(SrtDevError(['Conversion from vertex to facet',...
                    ' failed!\nThis is most likely due to numerical issues ',...
                    'in convhulln/CDDMEX in MPT.']));
                end        
            end
    end
    if options.verbose >= 2
        fprintf('Time to get inner min. H-rep  : %1.3f s\n',toc(timerVal));
    end

    % Store the facet version
    one_step_back_reach_polytope_underapprox = Polyhedron('H',...
        [one_step_back_reach_polytope_underapprox_A,...
         one_step_back_reach_polytope_underapprox_b], 'V',...
         one_step_back_reach_polytope_underapprox_V);

    if options.verbose == 3
        figure();
        if x_u_reaches_target_set.Dim <= 3
            before_projection_vrep_underapprox = ...
                Polyhedron('V',boundary_point_mat');            
            plot(x_u_reaches_target_set,'alpha',0.3);
            hold on;
            plot(before_projection_vrep_underapprox);
        else
            one_step_back_reach_polytope_underapprox_2D_v = ...
                Polyhedron('V',boundary_point_mat(1:2,:)');            
            hold on;
            plot(one_step_back_reach_polytope_underapprox_2D_v,'color','r',...
                'alpha',0.6);
            view([0 90]);
            box on;
        end
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% LRS-based one-step back reach set computation via Blachini/Nyugen's steps
function one_step_back_reach_polytope = lrsOneStepBackReachSet( ...
    current_time, sys, target_set, current_safe_set, verbose)

    if current_safe_set.isEmptySet() || target_set.isEmptySet()
        throwAsCaller(SrtInvalidArgsError('Recursion led to an empty set.'));
    end
    if ~target_set.isFullDim()
        disp(target_set)
        throwAsCaller(SrtInvalidArgsError('Cannot handle low-dim target set.'));
    end

    %% Compute the polytope in X*U 
    % Compute the polytope in X*U such that there is some x and u in the current
    % time which on application of dynamics will lie in target_set. We also add
    % the constraints defined by the current safe set as well as the input
    % space. This ensures that when we project 1) we will lie in a subset of
    % current_safe_set (effect of intersection after projection), and 2) the
    % input "selected" during the projection remains feasible.
    if sys.islti
        x_u_reaches_target_set_A = [target_set.A * [sys.state_mat sys.input_mat];
            [zeros(size(sys.input_space.A,1), sys.state_dim), sys.input_space.A];
            [current_safe_set.A zeros(size(current_safe_set.A,1), sys.input_dim)]];
        x_u_reaches_target_set_b = [target_set.b;
                                    sys.input_space.b;
                                    current_safe_set.b];
    else
        x_u_reaches_target_set_A = [target_set.A * ...
                [sys.state_mat(current_time) sys.input_mat(current_time)];
            [zeros(size(sys.input_space.A,1), sys.state_dim), sys.input_space.A];
            [current_safe_set.A zeros(size(current_safe_set.A,1), sys.input_dim)]];
        x_u_reaches_target_set_b = [target_set.b;
                                    sys.input_space.b;
                                    current_safe_set.b];
    end
    %% Compute projection
    if verbose >= 2
        timerVal = tic;
    end
 %    % Use GeoCalcLib's projection directly
 %    [one_step_back_reach_polytope_A_red,one_step_back_reach_polytope_b_red] =...
 %        projectPolyhedron(x_u_reaches_target_set_A, x_u_reaches_target_set_b,...
 %            uint32(sys.input_dim));
    one_step_back_reach_polytope_V_red = vertexEnumeration(...
        x_u_reaches_target_set_A, x_u_reaches_target_set_b);    
    if isempty(one_step_back_reach_polytope_V_red)
        throw(SrtRuntimeError('Recursion led to an empty target set!'));
    end
    % Project and reduce vertices
    one_step_back_reach_polytope_V = vertexReduction(...
        one_step_back_reach_polytope_V_red(:,1:sys.state_dim),...
        ones(size(one_step_back_reach_polytope_V_red,1),1));
    if verbose >= 2
        fprintf('Time to get v-rep             : %1.3f s\n',toc(timerVal));
        fprintf('Number of vertices            : %d vertices\n', ...
            length(one_step_back_reach_polytope_V_red));
    end    
    if verbose >= 2
        timerVal = tic;
    end
    [one_step_back_reach_polytope_A_red,one_step_back_reach_polytope_b_red] =...
        facetEnumeration(one_step_back_reach_polytope_V,...
        ones(size(one_step_back_reach_polytope_V,1),1));
    if verbose >= 2
        fprintf('Time to get H-rep             : %1.3f s\n',toc(timerVal));
        fprintf('Number of halfspaces          : %d inequalities\n', ...
            length(one_step_back_reach_polytope_b_red));
    end
    %% Compute min HRep()
    if verbose >= 2
        timerVal = tic;
    end
    % LRS based inequality reduction
    [one_step_back_reach_polytope_A,one_step_back_reach_polytope_b] = ...
        inequalityReduction(one_step_back_reach_polytope_A_red,...
            one_step_back_reach_polytope_b_red);
    one_step_back_reach_polytope = Polyhedron('H',...
        [one_step_back_reach_polytope_A, one_step_back_reach_polytope_b]);
    if verbose >= 2
        fprintf('Time to get minimal H-rep     : %1.3f s\n',toc(timerVal));
        fprintf('Number of minimal halfspaces  : %d inequalities\n', ...
                length(one_step_back_reach_polytope_b));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     %% Sanity check to make sure that vrep and hrep are consistent
%     p_v = Polyhedron('V',one_step_back_reach_polytope_underapprox_V);
%     p_h = Polyhedron('H',[one_step_back_reach_polytope_underapprox_A,...
%         one_step_back_reach_polytope_underapprox_b]);
%     p_h = Polyhedron('H',one_step_back_reach_polytope_underapprox.H);        
%     p_h.minHRep();
%     if verbose >= 2
%         timerVal = tic;
%     end
%     % Vertex representation contains half-space representation
%     flag1 = p_v.contains(p_h);   % Converts things to half-space btw
%     % Half-space representation contains vertex representation
%     flag2 = all(p_h.contains(p_v.V'));
%     if flag1 && flag2 
%         % all ok since p_v.eq(p_h) = 1
%         if verbose >= 2
%             disp('Both representations are consistent');
%         end
%     elseif flag1 || flag2 
%         warning('SReachTools:runTime',sprintf(['Strict equality not ',...  
%             'observed between the halfspace and vertex representation ',...
%             'of the polytope at the projection step for lag-under.',...
%             '\nHowever, we continue since one of the representations ',...
%             'contain the other.\n',...
%             'V-rep contains H-rep: %d\n',...
%             'H-rep contains V-rep: %d\n',...
%             'Ideally, this should not have happened! Please treat the ',...
%             'results with suspicions! This phenomena is most likely ',...
%             'due to numerical issues in convhulln/CDDMEX in MPT.'],...
%             flag1, flag2));
%     else
%         throw(MException('SReachTools:badConv','no-contain'));
%     end
%     if verbose >= 2
%         fprintf('Time for sanity check         : %1.3f s\n',toc(timerVal));
%     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Search for LRS to uncomment the LRS code.
