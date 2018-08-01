classdef TargetTube
% Create a target tube object
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
% Note that both of the above mechanisms will yield a target tube of length
% time_horizon+1 --- T_0, T_1, ..., T_{time_horizon}.
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
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
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
    properties
        dim
    end

    methods
        function obj = TargetTube(varargin)
        % TargetTube constructor
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
        %        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
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
                        exc = SrtInvalidArgsError(['Invalid arguments, ', ...
                            'see help.']);
                        throwAsCaller(exc);
                end

                validateattributes(target_set, {'Polyhedron'}, {'nonempty'});
                validateattributes(safe_set, {'Polyhedron'}, {'nonempty'});
                validateattributes(time_horizon, {'numeric'}, ...
                    {'integer', '>', 0});

                if safe_set.Dim ~= target_set.Dim
                    throwAsCaller(SrtInvalidArgsError(['Safe and target ', ...
                        'sets must be of the same dimension']));
                end

                obj.tube(1, time_horizon+1) = target_set;
                obj.dim = safe_set.Dim;
                for itt = 1:time_horizon
                    obj.tube(itt) = safe_set;
                end
            else
                validateattributes(varargin{nargin}, {'Polyhedron'}, ...
                        {'nonempty'})

                obj.tube(1, nargin) = varargin{nargin};
                obj.dim = obj.tube(nargin).Dim;
                for itt = 1:nargin-1
                    validateattributes(varargin{itt}, {'Polyhedron'}, ...
                        {'nonempty'})
                    if varargin{itt}.Dim ~= obj.dim
                        throwAsCaller(SrtInvalidArgsError(['All sets in the target tube must be of the same ', ...
                         'dimension']));
                    end

                    obj.tube(itt) = varargin{itt};
                end
            end
        end

        function varargout = subsref(obj, s)
        % Overloading of subsref
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
        % Inputs:
        % -------
        %   s - Subscripted reference struct array
        % 
        % see also subsref
        % 
        % Outputs: None
        % 
        % ====================================================================
        %
        %   This function is part of the Stochastic Optimal Control Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 

            if strcmp(s(1).type, '()')
                [varargout{1:nargout}] = subsref(obj.tube, s);
            else
                [varargout{1:nargout}] = builtin('subsref', obj, s); 
            end
        end

        function n_tubes = length(obj)
        % Overloading of length
        % ====================================================================
        %
        % Function to get length of target tube, returns length of tube
        % property
        %
        % Usage: Overloading of MATLAB length function
        % 
        % ====================================================================
        %
        % length(obj)
        % 
        % Inputs:  None
        % Outputs: None
        % 
        % ====================================================================
        %
        %   This function is part of the Stochastic Optimal Control Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 

            n_tubes = length(obj.tube);
        end

        function disp(obj)
        % Overloading of disp
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
        % Inputs:  None
        % Outputs: None
        % 
        % ====================================================================
        %
        %   This function is part of the Stochastic Optimal Control Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 

            fprintf('TargetTube of %d sets\n', length(obj));
        end

        function [concat_target_tube_A, concat_target_tube_b] = concat(obj)
        %  Get concatenated target tube
        % ======================================================================
        %
        % This method computes the concatenated target tube, 
        % safe_set^{time_horizon -1 } x target_set, a huge polyhedron in the
        % (obj.dim x time_horizon)-dimensional Euclidean space.
        % The output matrices satisfy the relation that the a concatenated 
        % state vector X lies in the reach-avoid tube if and only if
        % 
        % concat_target_tube_A * X <= concat_target_tube_b 
        %
        % Usage: See getFtLowerBoundTargetTube.
        %
        % ======================================================================
        %
        % [concat_target_tube_A, concat_target_tube_b] = concat(obj);
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
        % ======================================================================
        % 
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
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
        
        function [contains_flag] = contains(obj,X)
        %  Check if a given concatenates state
        % trajectory lies in the target tube
        % ======================================================================
        %
        % This method is a wrapper over MPT's Polyhedron/contains 
        % 
        % Usage: See checkViaMonteCarloSims
        %
        % ======================================================================
        %
        % [contains_flag] = contains(obj,X);
        % 
        % Inputs: 
        % -------
        %
        %   X             - Concatenated state vector (or collection of it
        %                   arranged columnwise)
        %    
        % Outputs:
        % --------
        %   contains_flag - Row vector of length equal to columns in X
        %
        % Notes:
        % ------
        % * This function is useful for Monte-Carlo simulation-based
        %   verification
        % 
        % ======================================================================
        % 
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        %
        %  

            if size(X, 1) ~= length(obj) * obj.dim
                throwAsCaller(SrtInvalidArgsError(['Concatenated state vector/matrix provided has incorrect ',...
                    'number of rows']));
            end

            if size(X, 2) <= 0
                throwAsCaller(SrtInvalidArgsError(['Concatenated state vector/matrix provided should have at ',...
                 'least one column']));
            end
            
            % Matrix to store the result for each time instant, trajectory
            contains_flag_all = zeros(length(obj),size(X,2));
            % Iterate over all time indexes since MPT contains will do
            % unnecessarily n_poly x n_state comparisons
            for t_indx = 1:length(obj)
                X_at_t_indx = X((t_indx-1)*obj.dim + 1: t_indx*obj.dim,:);
                contains_flag_all(t_indx,:) = obj.tube(t_indx).contains(X_at_t_indx);
            end
            contains_flag = min(contains_flag_all);
        end
        
        function v = end(obj, ~, ~)
        %  Overloading of MATLAB internal end 
        % function for objects
        % ======================================================================
        %
        % Overloading of MATLAB internal end function for objects. Allows
        % for use of end in array indexing
        % 
        % ======================================================================
        %
        % v = end(obj, ~, ~);
        % 
        % Inputs: None
        % Outputs:
        % --------
        %   v - Value of end; equal to the length of the target tube
        %
        % ======================================================================
        % 
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        %
        %
          
            v = length(obj);
        end
    end    
end
