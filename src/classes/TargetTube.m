classdef TargetTube
% SReachTools/TargetTube: Create a target tube object
% ==========================================================================
%
% Target tube class
%
% Usage:
% ------
% % Three different calling mechanisms, reach-avoid problem
% tt = TargetTube('reach-avoid', safe_set, target_set, time_horizon);
% 
% % Viability Problem
% tt = TargetTube('viability', safe_set, time_horizon);
% 
% % General tube
% % Can use general Polyhedron objects, for self-containment of the 
% % usage example will use empty Polyhedron objects
% tt = TargetTube(Polyhedron(), Polyhedron(), Polyhedron());
%   
% ==========================================================================
%
% TargetTube Properties:
% ------------------------
%
% TargetTube Methods:
% ---------------------
%   TargetTube/TargetTube - Class constructor
% 
% Notes:
% ------
% * MATLAB DEPENDENCY: MPT 3.0
% 
% =========================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/abyvinod/SReachTools/blob/master/LICENSE
% 
% 
    properties (Access = private)
        % TargetTube/tube
        % ==================================================================
        % 
        % Polyhedron array of sets for the target tube
        %
        tube = Polyhedron();
    end

    methods
        function obj = TargetTube(varargin)
        % SReachTools/TargetTube/TargetTube  TargetTube constructor
        % ====================================================================
        %
        % TargetTube class constructor
        %
        % Usage:
        % ------
        % % Three different calling mechanisms, reach-avoid problem
        % tt = TargetTube('reach-avoid', safe_set, target_set, time_horizon);
        % 
        % % Viability Problem
        % tt = TargetTube('viability', safe_set, time_horizon);
        % 
        % % General tube
        % % Can use general Polyhedron objects, for self-containment of the 
        % % usage example will use empty Polyhedron objects
        % tt = TargetTube(Polyhedron(), Polyhedron(), Polyhedron());
        % 
        % ====================================================================
        %
        % obj = TargetTube('reach-avoid', safe_set, target_set, time_horizon)
        % obj = TargetTube('viability', safe_set, time_horizon)
        % obj = TargetTube(poly1, ..., polyN)
        % 
        % Inputs:
        % -------
        %   safe_set          - Polyhedron object for set of safe states
        %   target_set        - Polyhedron object for set of target states
        %   time_hoizon       - Time horizon of target tube/reach-avoid problem
        %   poly1, ..., polyN - Polyhedron objects
        % 
        % Outputs:
        % --------
        %   obj - TargetTube object
        % 
        % ====================================================================
        %
        %   This function is part of the Stochastic Optimal Control Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/abyvinod/SReachTools/blob/master/LICENSE
        %   

            if ischar(varargin{1})
                % first argument character, specifying reach-avoid or viability
                switch(lower(varargin{1}))
                    case 'reach-avoid'
                        safe_set     = varargin{2};
                        target_set   = varargin{3};
                        time_horizon = varargin{4};
                    case 'viability'
                        safe_set     = varargin{2};
                        time_horizon = varargin{3};
                        target_set   = safe_set;
                    otherwise
                        error('Invalid arguments, see help.')
                end

                validateattributes(target_set, {'Polyhedron'}, {'nonempty'});
                validateattributes(safe_set, {'Polyhedron'}, {'nonempty'});
                validateattributes(time_horizon, {'numeric'}, ...
                    {'integer', '>', 0});

                assert(safe_set.Dim == target_set.Dim,...
                   'SReachTools:invalidArgs',...
                   'Safe and target sets must be of the same dimension');

                obj.tube(1, time_horizon) = target_set;
                for itt = 1:time_horizon-1
                    obj.tube(itt) = safe_set;
                end
            else
                validateattributes(varargin{nargin}, {'Polyhedron'}, ...
                        {'nonempty'})

                obj.tube(1, nargin) = varargin{nargin};
                set_dim = obj.tube(nargin).Dim;
                for itt = 1:nargin-1
                    validateattributes(varargin{itt}, {'Polyhedron'}, ...
                        {'nonempty'})
                    assert(varargin{itt}.Dim == set_dim, ...
                        'SReachTools:invalidArgs', ...
                        ['All sets in the target tube must be of the same ', ...
                         'dimension']);

                    obj.tube(itt) = varargin{itt};
                end
            end
        end

        function varargout = subsref(obj, s)
        % SReachTools/TargetTube/subsref  Overloading of subsref
        % ====================================================================
        %
        % Overloading of builting MATLAB subsref for the class. If the first
        % subscripted reference for the object is of type '()' then the 
        % remaining references are made on the tube (Polyhedron object array) 
        % property, otherwise the TargetTube object is refernced.
        %
        % Usage: Overloading of MATLAB internal
        % 
        % ====================================================================
        %
        % subsref(obj, s)
        % 
        % Input:
        % ------
        %   s - Subscripted reference struct array
        % 
        % see also subsref
        % 
        % Outpus: None
        % 
        % ====================================================================
        %
        %   This function is part of the Stochastic Optimal Control Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/abyvinod/SReachTools/blob/master/LICENSE
        % 

            if strcmp(s(1).type, '()')
                [varargout{1:nargout}] = subsref(obj.tube, s);
            else
                [varargout{1:nargout}] = builtin('subsref', obj, s); 
            end
        end

        function n_tubes = length(obj)
        % SReachTools/TargetTube/length  Overloading of length
        % ====================================================================
        %
        % Function to get length of target tube, returns length of tube
        % property
        %
        % Usage: Overloadig of MATLAB length function
        % 
        % ====================================================================
        %
        % length(obj)
        % 
        % Input:  None
        % Outpus: None
        % 
        % ====================================================================
        %
        %   This function is part of the Stochastic Optimal Control Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/abyvinod/SReachTools/blob/master/LICENSE
        % 

            n_tubes = length(obj.tube);
        end

        function disp(obj)
        % SReachTools/TargetTube/disp  Overloading of disp
        % ====================================================================
        %
        % Overloading of builting MATLAB disp for the class.
        %
        % Usage: Overloading of MATLAB internal
        % 
        % ====================================================================
        %
        % disp(obj)
        % 
        % Input:  None
        % Outpus: None
        % 
        % ====================================================================
        %
        %   This function is part of the Stochastic Optimal Control Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/abyvinod/SReachTools/blob/master/LICENSE
        % 

            disp(sprintf('   TargetTube of %d sets\n', length(obj)));
        end

        function [concat_target_tube_A, concat_target_tube_b] = concat(obj)
        % SReachTools/TargetTube/concat: Get concatenated target tube
        % ============================================================================
        %
        % This method computes the concatenated target tube, 
        % safe_set^{time_horizon -1 } x target_set, a huge polyhedron in the
        % (sys.state_dimension x time_horizon)-dimensional Euclidean space.
        % The output matrices satisfy the relation that the a concatenated 
        % state vector X lies in the reach-avoid tube if and only if
        % 
        % concat_target_tube_A * X <= concat_target_tube_b 
        %
        % Usage: See getFtLowerBoundTargetTube.
        %
        % ============================================================================
        %
        % [concat_target_tube_A, concat_target_tube_b] = concaat(obj);
        % 
        % Inputs: None
        %    
        % Outputs:
        % --------
        %   concat_target_tube_A - State matrix concatenated for target tube
        %   concat_target_tube_b - Input matrix concatenated for target tube
        %
        % Notes:
        % ------
        % * This function also serves as a delegatee for input handling.
        % 
        % ============================================================================
        % 
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://github.com/abyvinod/SReachTools/blob/master/LICENSE
        %
        %

            %% Construction of the concatenated target tube
            tube_A_mats = cell(1, length(obj));
            [tube_A_mats{:}] = obj.tube(:).A;
            concat_target_tube_A = blkdiag(tube_A_mats{:});

            tube_b_vecs = cell(1, length(obj));
            [tube_b_vecs{:}] = obj.tube(:).b;
            concat_target_tube_b = vertcat(tube_b_vecs{:});
        end
    end
end