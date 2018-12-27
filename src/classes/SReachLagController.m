classdef SReachLagController
% A feasible controller that maximizes stochastic reachability via 
% Lagrangian-based computations    
    properties
        % SReachLagController/tube
        % ==================================================================
        % 
        % Polyhedron array of sets for the tube (half-space representation
        % preferred)
        %
        tube = Polyhedron();
        % SReachLagController/system
        % ==================================================================
        % 
        % LtvSystem or LtiSystem object describing the dynamics
        %
        system = [];
        % SReachLagController/dist_set
        % ==================================================================
        % 
        % Polyhedron or SReachEllipsoid object describing the disturbance set
        %
        dist_set = [];
        % SReachLagController/time_horizon
        % ==================================================================
        % 
        % Time horizon up to which the controller has been defined 
        % (Defined using the length(obj.tube) and must be >0)
        %
        time_horizon
    end

    methods
        function obj = SReachLagController(sys, dist_set, tube)
            % Add the system after input handling
            validateattributes(sys, {'LtiSystem','LtvSystem'}, {'nonempty'});
            if sys.input_dim < 1
                 throwAsCaller(SrtInvalidArgsError(['Requires',...
                     ' an uncontroller system.']));
            end
            obj.system = sys;
            % Add the disturbance after input handling (Can be empty if
            % unperturbed)
            validateattributes(dist_set, {'Polyhedron','SReachEllipsoid'},...
                {});
            if sys.dist_dim > 0
                switch class(dist_set)
                    case 'Polyhedron'
                        if dist_set.Dim ~= sys.dist_dim
                            throwAsCaller(SrtInvalidArgsError(['Mismatch in',...
                                ' dimension of the disturbance set and the ',...
                                'additive noise in the system.']));
                        end
                        if dist_set.isEmptySet()
                            throwAsCaller(SrtInvalidArgsError(['Requires a',...
                                ' non-empty disturbance set.']));
                        end
                    case 'SReachEllipsoid'
                        if dist_set.dim ~= sys.dist_dim
                            throwAsCaller(SrtInvalidArgsError(['Mismatch in',...
                                ' dimension of the disturbance set and the ',...
                                'additive noise in the system.']));
                        end
                end
                obj.dist_set = dist_set;
            else
                obj.dist_set = [];
            end
            % Add the target tube after input handling            
            validateattributes(tube, {'Polyhedron','Tube'}, {'nonempty'});
            tube = Tube.polyArray2Tube(tube);
            if tube.dim ~= sys.state_dim
                throwAsCaller(SrtInvalidArgsError(['Mismatch in ',...
                    'dimension of the collection of time-stamped stochastic',...
                    ' reach set underapproximations and the state.']));
            end
            obj.tube = tube;
            % Time horizon
            obj.time_horizon = length(tube) - 1;
            if obj.time_horizon < 1
                throwAsCaller(SrtInvalidArgsError('Time horizon must be > 0.'));
            end
        end
        
        function action = getInput(obj, current_state, current_time)
            validateattributes(current_time, {'numeric'}, {'scalar', ...
                'integer', '>=', 0, '<=', obj.time_horizon - 1});
            validateattributes(current_state, {'numeric'}, {'column', ...
                'size', [obj.system.state_dim 1]});
            % time goes from 0 to N, and tube_indx goes from 1 to N+1
            tube_indx = current_time + 1;
            
            % Return NaN if invalid current state
            if ~(obj.tube(tube_indx).contains(current_state)) ||...
                    any(isnan(current_state))
                throwAsCaller(SrtInvalidArgsError('Invalid current state'));
%                 action = nan(obj.system.input_dim, 1);
%                 return
            end
            
            % Ensure that next effective target set has a half-space
            % representation
            target_set = obj.tube(tube_indx+1);
            if ~target_set.hasHRep
                target_set.computeHRep();
            end
            
            % Compute the effective input set
            %       (R_{t+1} \ominus F_t*E_t) \oplus {-Ax_t}
            if ~isempty(obj.dist_set)
                scaled_dist_set = obj.system.dist_mat(current_time) * obj.dist_set;

                effective_input_set = Polyhedron('H', ...
                    [target_set.A * obj.system.input_mat(current_time),...
                        (target_set.b -...
                            scaled_dist_set.support(target_set.A') -...
                            target_set.A * ...
                            obj.system.state_mat(current_time) * current_state);
                     obj.system.input_space.A, obj.system.input_space.b]);
            else
                effective_input_set = Polyhedron('H', ...
                    [target_set.A, target_set.b -...
                        target_set.A*obj.system.state_mat(current_time)*current_state;
                     obj.system.input_space.A, obj.system.input_space.b]);
            end
            
            % Compute a feasible action via MPT's interior point (Chebyshev
            % centering)
            if effective_input_set.isEmptySet()
                action = nan(obj.system.input_dim, 1);
            else
                sol = effective_input_set.interiorPoint();
                action = sol.x;
            end
        end
    end    
end
