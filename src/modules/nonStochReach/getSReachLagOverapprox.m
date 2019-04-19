function varargout = getSReachLagOverapprox(sys, target_tube,...
    disturbance_set, options)
% Get the overapproximation of the stochastic reach set
% ============================================================================
%
% This function will compute the overapproximation of the stochastic reach
% set via Algorithm 2 in
% 
%      J. D. Gleason, A. P. Vinod, and M. M. K. Oishi. 2018. Lagrangian 
%      Approximations for Stochastic Reachability of a Target Tube. 
%      online. (2018). https://arxiv.org/abs/1810.07118
%
% Usage: See examples/lagrangianApproximations.m
%   
% ============================================================================
%
% [overapprox_set, overapprox_tube] = getSReachLagUnderapprox(sys,...
%       target_tube, disturbance_set)
%
% Inputs:
% -------
%   sys             - LtiSystem object
%   target_tube     - Tube object 
%   disturbance_set - Polyhedron or SReachEllipsoid object (bounded disturbance 
%                     set)
%   options         - Struct of reach set options, see SReachSetOptions
%
% Outputs:
% --------
%   overapprox_set - Polyhedron object for the overapproximation of the 
%                    stochastic reach set
%   overapprox_tube- [Available for 'VHmethod' only] Tube comprising of an
%                    overapproximation of the stochastic reach sets across the
%                    time horizon
%
% Notes:
% * From polyhedral computation theory (implemented here when
%   options.compute_style is VHmethod), intersections and Minkowski differences
%   are best performed in facet representation and Minkowski sums are best
%   performed in vertex representation. However, since in this computation, all
%   three operations are required, scalability of the algorithm is severly
%   hampered, despite theoretical elegance.
% * Using support functions, an arbitrarily tight polytopic overapproximation of
%   the set may be computed via convex optimization (linear or second order-cone
%   programs).
%
% ============================================================================
%
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://sreachtools.github.io/license/
%

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
    
    % Check if prob_str and method_str of options are consistent        
    if ~strcmpi(options.prob_str, 'term')
        throwAsCaller(...
            SrtInvalidArgsError('Mismatch in prob_str in the options'));
    end
    if ~strcmpi(options.method_str, 'lag-over')
        throwAsCaller(...
            SrtInvalidArgsError('Mismatch in method_str in the options'));
    end            
    
    switch lower(options.compute_style)
        case 'vfmethod'
            %% Use vertex-facet enumeration requiring recursion
            [effective_target_tube] = computeViaRecursion(sys, target_tube,...
                disturbance_set, options);
            varargout{1} = effective_target_tube(1);
            varargout{2} = effective_target_tube;
        case 'support'
            %% Use recursion-free support function approach
            [effective_target_set] = computeViaSupportFn(sys, target_tube,...
                disturbance_set, options);
            varargout{1} = effective_target_set;            
        otherwise
            throwAsCaller(SrtInvalidArgsError(['Invalid computation style ', ...
                'specified']));
    end
end

function [effective_target_tube] = computeViaRecursion(sys, target_tube,...
    disturbance_set, options)   
% This private function implements the recursion-based computation which
% internally requires vertex-facet enumeration => performs really well for
% low dimension systems but scales very poorly
%
% ============================================================================
%
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://sreachtools.github.io/license/
%

    tube_length = length(target_tube);
    if sys.islti()
        inverted_state_matrix = inv(sys.state_mat);
        minus_bu = (-sys.input_mat) * sys.input_space;
        minus_scaled_dist_set = (-sys.dist_mat) * disturbance_set;
    end
    if options.verbose >= 1
        fprintf('Time_horizon: %d\n', tube_length-1);
    end

    effective_target_tube = repmat(Polyhedron(), tube_length, 1);
    effective_target_tube(end) = target_tube(end);
    if tube_length > 1
        for itt = tube_length-1:-1:1
            current_time = itt - 1;
            if options.verbose >= 1
                fprintf('Computation for time step: %d\n', current_time);
            end
            if sys.isltv()
                % Overwrite the following parameters with their
                % time-varying counterparts
                inverted_state_matrix = inv(sys.state_mat(current_time));
                minus_bu = (-sys.input_mat(current_time)) * sys.input_space;
                minus_scaled_dist_set = (-sys.dist_mat(current_time)) *...
                    disturbance_set;
            end
            
            if isa(disturbance_set,'Polyhedron') && disturbance_set.isEmptySet
                % No augmentation
                new_target = effective_target_tube(itt+1);
            elseif ~effective_target_tube(itt+1).isEmptySet()
                % Compute a new target set for this iteration that is robust to 
                % the disturbance
                new_target = minus_scaled_dist_set.plus(...
                    effective_target_tube(itt+1));
            else
                throw(SrtRuntimeError('Recursion led to an empty target set!'));
            end

            % One-step backward reach set
            one_step_backward_reach_set = inverted_state_matrix * ...
                (new_target + minus_bu);

            % Guarantee staying within target_tube by intersection
            effective_target_tube(itt) = intersect(...
                one_step_backward_reach_set, target_tube(itt));
        end
    end
end

function [effective_target_set] = computeViaSupportFn(sys, target_tube,...
    disturbance_set, options)   
% This private function implements the support function-based
% recursion-free implementation of Lagrangian overapproximation of the
% stochastic reach set.
%
% ============================================================================
%
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://sreachtools.github.io/license/
%

    if isempty(options.equi_dir_vecs)
        throwAsCaller(SrtInvalidArgsError(['Expected non-empty ',...
            'equi_dir_vecs. Faulty options structure provided!']));
    end
    % Get size of equi_dir_vecs
    [dir_vecs_dim, n_vertices] = size(options.equi_dir_vecs);

    if dir_vecs_dim ~= (sys.state_dim) || n_vertices < 3
        throwAsCaller(SrtInvalidArgsError(['Expected (sys.state_dim + ',...
            'sys.input_dim)-dimensional collection of column vectors. ',...
            'Faulty options structure provided!']));        
    end

    effective_target_set_A = zeros(n_vertices, dir_vecs_dim);
    effective_target_set_b = zeros(n_vertices, 1);
    if options.verbose == 1
        fprintf('Evaluating support function: 00000/%5d',n_vertices);
    end
        
    for dir_indx = 1:n_vertices
        if options.verbose == 1
            fprintf('\b\b\b\b\b\b\b\b\b\b\b%5d/%5d',dir_indx, n_vertices);
        elseif options.verbose == 2
            fprintf('\n\nEvaluating support function: %5d/%5d\n\n',...
                dir_indx, n_vertices);
        end
        ell = options.equi_dir_vecs(:, dir_indx);
        effective_target_set_A(dir_indx, :)= ell';
        effective_target_set_b(dir_indx) =...
            support(ell, sys, target_tube, disturbance_set, options);
        if options.verbose >= 3 && size(ell,1) <=3
            figure(200);
            hold on;
            title(sprintf('support %d/%d',dir_indx, n_vertices));
            p_temp = Polyhedron('H',[ell' effective_target_set_b(dir_indx)]);
            plot(p_temp.intersect(target_tube(1)),'alpha',0);
            %scatter(ell(1),ell(2),300,'rx');
            quiver(0,0,ell(1),ell(2),'rx');
            axis equal;
            drawnow;
        end
    end
    switch options.vf_enum_method
        case 'cdd'                                      % MPT/CDD approach
            effective_target_set = Polyhedron('H', [effective_target_set_A, ...
                effective_target_set_b; target_tube(1).A, target_tube(1).b]);
            effective_target_set.minHRep();
        case 'lrs'                                      % LRS approach
            effective_target_set_A_red = [effective_target_set_A; 
                target_tube(1).A];
            effective_target_set_b_red = [effective_target_set_b; 
                target_tube(1).b];
            [effective_target_set_A_irred, effective_target_set_b_irred] = ...
                inequalityReduction(effective_target_set_A_red, ...
                    effective_target_set_b_red);
            effective_target_set = Polyhedron('H', ...
                [effective_target_set_A_irred, effective_target_set_b_irred]);
    end
    if options.verbose >= 1
        fprintf('\n');
    end
end

function [val] = support(ell, sys, target_tube, dist_set, options)
% This private function implements the support function-based
% recursion-free implementation of Lagrangian overapproximation of the
% stochastic reach set. Specifically, it is called by computeViaSupportFn for
% each sample of the support function.
%
% ============================================================================
%
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://sreachtools.github.io/license/
%

    % Time horizon
    time_horizon = length(target_tube) - 1;
    % No. of facets in the input space
    n_lin_input = size(sys.input_space.A,1);
    % Half-space representation for the concatenated target tube
    [concat_target_tube_A, concat_target_tube_b] = target_tube.concat();

    cvx_begin
        if options.verbose >= 2
            cvx_quiet false
        else
            cvx_quiet true
        end
            
        % Dummy variable (no. of time steps-long)
        variable nu(sys.state_dim, time_horizon);
        % Dual variable for the target set
        variable dual_var_target(size(concat_target_tube_A,1), 1) nonnegative;
        % Dual variable for the input space
        variable dual_var_input(n_lin_input, time_horizon) nonnegative;        
        

        switch class(dist_set)
            case 'Polyhedron'
                % Dual variable for the disturbance set
                variable dual_var_dist(size(dist_set.A,1),time_horizon) nonnegative;  
                % Since the input set and disturbance set are time-invariant,
                % we can compute b' * [z_0 ... z_{N-1}] and then sum it up!
                minimize ((concat_target_tube_b'*dual_var_target) + sum(sys.input_space.b' * dual_var_input) + sum(dist_set.b'* dual_var_dist));
            case 'SReachEllipsoid'
                % Slack variable to account for Minkowski sum of input space and
                % disturbance set
                variable slack_var_inputdist(time_horizon,1);
                % Here, we have a closed form expression for the support
                % function of the disturbance set => Use slack variables to
                % enforce the constraint in an epigraph form
                minimize ((concat_target_tube_b'*dual_var_target) + sum(slack_var_inputdist))
            otherwise
                throwAsCaller(SrtInvalidArgsError(sprintf(['Disturbance (%s',...
                    ') is not configured as of yet'], class(scaled_dist_set))));
        end
        
        
        subject to            
            dual_var_start_indx = 1;
            dual_var_end_indx = size(target_tube(1).A,1);
                
            for tube_indx = 1:time_horizon + 1
                % current_time, denoted by k, is t_indx-1 and goes from 0 to N
                current_time = tube_indx - 1;
                % Computation of A_k^{-1}
                inv_sys_now = inv(sys.state_mat(current_time));
                if current_time >= 1
                    % Computation of A_{k-1}^{-1}
                    inv_sys_prev = inv(sys.state_mat(current_time - 1));                
                end

                %% Target set at t \in N_{[0, N]}
                if tube_indx == 1
                    % (z_Target_0)' * A_Target_0 == (l - nu_0)'
                    (dual_var_target(dual_var_start_indx:dual_var_end_indx))'... 
                        * (target_tube(1).A) == (ell - nu(:,1))';                    
                elseif tube_indx < time_horizon + 1
                    % (z_Target_k)' * A_Target_k ==
                    %                       (nu_{k-1}'*A_sys_(k-1)^{-1} - nu_k')
                    (dual_var_target(dual_var_start_indx:dual_var_end_indx))'...
                        * (target_tube(tube_indx).A) == ...
                            (nu(:,tube_indx - 1)' * inv_sys_prev - ...
                                nu(:, tube_indx)');
                elseif tube_indx == time_horizon + 1
                    % A_Target_{N}'*z_Target_{N} == A_sys_(N-1)^{-T}*nu_{N-1}
                    % N + 1 is the corresponding tube_indx for k=N
                    (dual_var_target(dual_var_start_indx:dual_var_end_indx))'...
                        * (target_tube(time_horizon + 1).A) == ...
                            nu(:,time_horizon)' * inv_sys_prev;
                end
                
                % Increment the start counter
                if tube_indx <= time_horizon
                    dual_var_start_indx = dual_var_end_indx + 1;
                    dual_var_end_indx = dual_var_start_indx - 1 + ...
                         size(target_tube(tube_indx+1).A,1);
                end
                
                %% Support function of (-BU) + (-FE) for k from 0 to N-1
                if tube_indx <= time_horizon
                    % z_u_k' * A_u == - nu_k' * A_k^{-1} * B_k
                    dual_var_input(:, tube_indx)' * sys.input_space.A ==...
                        - nu(:, tube_indx)' * inv_sys_now *...
                            sys.input_mat(current_time);
                    
                    switch class(dist_set)
                        case 'SReachEllipsoid'
                            % Compute F_k E for k from 0 to N-1
                            minus_fe_now = - sys.dist_mat(current_time) *...
                                dist_set;
                    
                            % z_u_k' * b_u + minus_fe_now.support(A^{-T}_k nu_k)
                            %                                  <= s_inputdist_k
                            dual_var_input(:, tube_indx)' * sys.input_space.b...
                                 + minus_fe_now.support(...
                                    inv_sys_now'*nu(:, tube_indx))...
                                    <= slack_var_inputdist(tube_indx);
                        case 'Polyhedron'
                            % A_E' * z_w_k == - nu_k' * A_k^{-1} * F_k
                            dual_var_dist(:,tube_indx)' * dist_set.A ==...
                                - nu(:, tube_indx)' * inv_sys_now * ...
                                    sys.dist_mat(current_time);
                        otherwise
                            throwAsCaller(SrtInvalidArgsError(sprintf( ...
                                'Unconfigured disturbance (%s) provided',...
                                class(scaled_dist_set))));
                    end
                end
            end
    cvx_end
    switch cvx_status
        case 'Solved'
            val = cvx_optval;
        case 'Inaccurate/Solved'
            val = cvx_optval;
            % TODO: Check if the solution is indeed accurate
        otherwise
            err_mesg = sprintf(['Support function computation via CVX ',...
                'failed. CVX status: %s'], cvx_status);
            throw(SrtRuntimeError(err_mesg));
    end
end
