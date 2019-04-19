classdef SReachLagController
% Controller that maximizes stochastic reachability via Lagrangian-based
% computations and respects hard control bounds
% ==========================================================================
%
% SReachLagController class
%
% Usage:
% ------
% % Given a system sys, probability threshold prob_thresh, and a target_tube, we
% % can compute the Lagrangian-based underapproximation via these two commands.
% lagunder_options = SReachSetOptions('term', 'lag-under', ...
%      'bound_set_method', 'ellipsoid', 'compute_style','vfmethod');
% [polytope_lagunder, extra_info_under] = SReachSet('term', 'lag-under', ...
%         sys, prob_thresh, target_tube, lagunder_options);
% % We compute the associated controller using the following command
% srlcontrol = SReachLagController(sys, bounded_dist_set, robust_reach_tube);
% 
% See also, example/cwhSReachSet.
%
% ==========================================================================
%
% SReachLagController Properties:
% -------------------------------
%   system          - System under study [LtvSystem/LtiSystem object]
%   dist_set        - Bounded disturbance set for which robustness computation
%                     has been completed [Polyhedron/SReachEllipsoid object]
%   tube            - Time-stamped robust sets in the state space from which
%                     there is a control that works irrespective of the
%                     disturbance in dist_set [Tube object]. See notes.
%   time_horizon    - Computed from the obj.tube (Scalar)
%
% SReachLagController Methods:
% ----------------------------
%   SReachLagController/SReachLagController 
%                   - Constructor for SReachLagController
%   SReachLagController/getInput 
%                   - Get the input at time t given the state at time t, that
%                     lies in the effective target tube
% 
% Notes:
% ------
% * MATLAB DEPENDENCY: MPT 3.0
% * Robust target tube can also be provided as an array of polyhedron.
%   The construct converts it into a tube object using
%   Tube.polyArray2Tube().
% 
% =========================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://sreachtools.github.io/license/
% 
% 
    properties
        % SReachLagController/tube
        % ==================================================================
        % 
        % Time-stamped robust sets in the state space from which there is a
        % control that works irrespective of the disturbance in dist_set [Tube
        % object]. When Polyhedron array of sets for the tube is given, the
        % constructor will use Tube.polyArray2Tube() to convert it into a tube.
        % Half-space representation preferred.
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
        % Constructor for SReachLagController 
        % ======================================================================
        % Inputs:
        % -------
        %   sys      - System under study [LtvSystem/LtiSystem object]
        %   dist_set - Bounded disturbance set for which robustness computation
        %              has been completed [Polyhedron/SReachEllipsoid object]
        %   tube     - Time-stamped robust sets in the state space from which
        %              there is a control that works irrespective of the
        %              disturbance in dist_set [Tube object]. See notes.
        %
        % Outputs:
        % --------
        %   obj      - SReachLagController object
        %
        % Notes:
        % ------
        % * Robust target tube can also be provided as an array of polyhedron.
        %   The construct converts it into a tube object using
        %   Tube.polyArray2Tube().
        %
        % ======================================================================
        % 
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://sreachtools.github.io/license/
        %
        %
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
            % Ensure that (t+1)-effective target set R_{t+1} has a half-space
            % representation
            for tube_indx = 1:length(tube)
                if ~tube(tube_indx).hasHRep
                    tube(tube_indx).minHRep();
                end
            end
            obj.tube = tube;
            
            % Time horizon
            obj.time_horizon = length(tube) - 1;
            if obj.time_horizon < 1
                throwAsCaller(SrtInvalidArgsError('Time horizon must be > 0.'));
            end
        end
        
        function action = getInput(obj, current_state, current_time)
        % Get the input at time t given the state at time t, that lies in the
        % effective target tube
        % ======================================================================
        % Inputs:
        % -------
        %   obj             - SReachLagController object
        %   current_state   - Current state (a obj.system.state_dim x 1 vector)
        %   current_time    - Current time (a scalar integer) 
        %
        % Outputs:
        % --------
        %   action          - Action to apply at time t
        %
        % Notes:
        % ------
        % * If infeasible initial state is provided or the feasible input space
        %   turns out to be empty, an Invalid Arguments error is thrown.
        %   Therefore, while using this function to generate trajectories, make
        %   sure that the disturbances lie in obj.dist_set.
        %
        % ======================================================================
        % 
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://sreachtools.github.io/license/
        %
        %
            validateattributes(current_state, {'numeric'}, {'column', ...
                'size', [obj.system.state_dim 1]});
            validateattributes(current_time, {'numeric'}, {'scalar', ...
                'integer', '>=', 0, '<=', obj.time_horizon - 1});
            % time goes from 0 to N, and tube_indx goes from 1 to N+1
            tube_indx = current_time + 1;
            if any(isnan(current_state)) ||...
                ~(obj.tube(tube_indx).contains(current_state))
                disp(current_state);
                disp(obj.tube(tube_indx));
                disp(obj.tube(tube_indx).contains(current_state));
                throwAsCaller(SrtInvalidArgsError('Invalid current state'));
            end
            
            % Target set is the next time-step 
            target_set = obj.tube(tube_indx + 1);            
            % Ensure that the effective target set has a half-space
            % representation
            if ~ target_set.hasHRep
                target_set.computeHRep();
            end

            % Feasible input set is (R_{t+1} \ominus F_t*E_t) \oplus {-A_t x_t}
            % where R_{t+1} is the (t+1)-effective target set, F_t is the
            % disturbance matrix, E_t is the disturbance set, A_t is the state
            % matrix, and x_t is the current state
            if ~isempty(obj.dist_set)
                % Compute F_t*E_t
                scaled_dist_set = obj.system.dist_mat(current_time) *...
                    obj.dist_set;

                % Compute the feasible set in the half space
                effective_input_set = Polyhedron('H', ...
                    [target_set.A * obj.system.input_mat(current_time),...
                        (target_set.b -...
                            scaled_dist_set.support(target_set.A') -...
                            target_set.A * ...
                            obj.system.state_mat(current_time) * current_state);
                     obj.system.input_space.A, obj.system.input_space.b]);
            else
                % Feasible input set is R_{t+1} \oplus {-A_t x_t}
                effective_input_set = Polyhedron('H', ...
                    [target_set.A, target_set.b - target_set.A * ...
                        obj.system.state_mat(current_time)*current_state;
                    obj.system.input_space.A, obj.system.input_space.b]);
            end
            
            % Compute a feasible action via MPT's interior point (Chebyshev
            % centering)
            if effective_input_set.isEmptySet()
                throw(SrtRuntimeError('Feasible input set is empty'));
            else
                sol = effective_input_set.interiorPoint();
                action = sol.x;
            end
        end
    end    
end
