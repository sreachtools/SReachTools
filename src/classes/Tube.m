classdef Tube
% Create a tube object
% ==========================================================================
%
% Tube class
%
% Usage:
% ------
% % Three different calling mechanisms, reach-avoid problem
% tt = Tube('reach-avoid', safe_set, target_set, time_horizon);
% 
% % Viability Problem
% tt = Tube('viability', safe_set, time_horizon);
%
% % Note that both of the above mechanisms will yield a tube of length
% % time_horizon+1 --- T_0, T_1, ..., T_{time_horizon}.
% 
% % General tube
% % Can use general Polyhedron objects, for self-containment of the 
% % usage example will use empty Polyhedron objects
% tt = Tube(Polyhedron(), Polyhedron(), Polyhedron());
%
% Given a tube 'orig_tube' and a polytope 'int_polytope', you can also 
% construct a tube that is the result of the intersection of the polytope 
% and sets in the tube as
% tt = Tube('intersect', orig_tube, int_polytope);
%   
% ==========================================================================
%
% Tube Properties:
% ------------------------
%   dim - Dimension of the tube, i.e. dimension of the sets in the tube
%
% Tube Methods:
% ---------------------
%   Tube/Tube       - Class constructor
%   concat          - Get concatenated tube (Cartesian product of the polytopes)
%   contains        - Check if a given concatenates state trajectory lies in the
%                     tube
%   polyArray2Tube  - Convert the array of polyhedra to a Tube object
% 
% Apart from these, the following MATLAB functions have been overloaded
%   disp            - Display the tube
%   length          - Provide the length of the tube 
%   end             - Indexing the end of the tube
%   subsref         - Permits use of indexing of tube
%
% Notes:
% ------
% * MATLAB DEPENDENCY: MPT 3.0
% 
% =========================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://sreachtools.github.io/license/
% 
% 
    properties (Access = private)
        % Tube/tube
        % ==================================================================
        % 
        % Polyhedron array of sets for the tube
        %
        tube = Polyhedron();
    end
    properties
        % Tube/dim
        % ==================================================================
        % 
        % Dimension of the tube / sets in the tube
        % 
        dim
    end

    methods
        function obj = Tube(varargin)
        % Tube constructor
        % ====================================================================
        %
        % Tube class constructor
        %
        % Usage:
        % ------
        % % Three different calling mechanisms, reach-avoid problem
        % tt = Tube('reach-avoid', safe_set, target_set, time_horizon);
        % 
        % % Viability Problem
        % tt = Tube('viability', safe_set, time_horizon);
        % 
        % % General tube
        % % Can use general Polyhedron objects, for self-containment of the 
        % % usage example will use empty Polyhedron objects
        % tt = Tube(Polyhedron(), Polyhedron(), Polyhedron());
        % 
        % Given a tube 'orig_tube' and a polytope 'int_polytope', you can also
        % construct a tube that is the result of the intersection of the
        % polytope and sets in the tube as
        % tt = Tube('intersect', orig_tube, int_polytope);
        %   
        % ====================================================================
        %
        % obj = Tube('reach-avoid', safe_set, target_set, time_horizon)
        % obj = Tube('viability', safe_set, time_horizon)
        % obj = Tube(poly1, ..., polyN)
        % obj = Tube('intersect', orig_tube, int_polytope);
        % 
        % Inputs:
        % -------
        %   safe_set          - Polyhedron object for set of safe states
        %   target_set        - Polyhedron object for set of target states
        %   time_hoizon       - Time horizon of tube/reach-avoid problem
        %   poly1, ..., polyN - Polyhedron objects
        %   orig_tube         - A tube object that needs to be intersected
        %   int_polytope      - Polyhedron object that is used for intersection
        % 
        % Outputs:
        % --------
        %   obj - Tube object
        % 
        % ====================================================================
        %
        %   This function is part of the Stochastic Optimal Control Toolbox.
        %   License for the use of this function is given in
        %        https://sreachtools.github.io/license/
        %   

            if isempty(varargin)
                err = SrtInvalidArgsError('Expected at least one input');
                throwAsCaller(err);
            end
            if ischar(varargin{1})
                % first argument character, specifying reach-avoid or viability
                tube_type = lower(varargin{1});
                switch(tube_type)
                    case 'reach-avoid'
                        safe_set     = varargin{2};
                        target_set   = varargin{3};
                        time_horizon = varargin{4};

                        validateattributes(target_set, {'Polyhedron'}, ...
                            {'nonempty'});
                        validateattributes(safe_set, {'Polyhedron'}, ...
                            {'nonempty'});
                        validateattributes(time_horizon, {'numeric'}, ...
                            {'integer', '>', 0});
                        
                        if safe_set.Dim ~= target_set.Dim
                            throwAsCaller(SrtInvalidArgsError(['Safe and', ...
                                ' target sets must be of the same dimension']));
                        end
                    case 'viability'
                        safe_set     = varargin{2};
                        time_horizon = varargin{3};
                        target_set   = safe_set;
                        
                        validateattributes(safe_set, {'Polyhedron'}, ...
                            {'nonempty'});
                        validateattributes(time_horizon, {'numeric'}, ...
                            {'integer', '>', 0});

                    case 'intersect'
                        orig_tube = varargin{2};
                        int_polytope = varargin{3};
                        
                        validateattributes(orig_tube, {'Tube'}, {'nonempty'});
                        validateattributes(int_polytope, {'Polyhedron'}, ...
                            {'nonempty'});                        
                        if orig_tube.dim ~= int_polytope.Dim
                            throw(SrtInvalidArgsError(['Dimension mismatch ', ...
                                'between intersecting polytope and tube.']));
                        end
                    otherwise
                        err = SrtInvalidArgsError(['Invalid arguments, ', ...
                            'see help.']);
                        throwAsCaller(err);
                end

                if strcmpi(tube_type,'viability') ||...
                    strcmpi(tube_type,'reach-avoid')
                    obj.tube(1, time_horizon+1) = target_set;
                    obj.dim = safe_set.Dim;
                    for itt = 1:time_horizon
                        obj.tube(itt) = safe_set;
                    end
                elseif strcmpi(tube_type,'intersect')                    
                    for itt = 1:length(orig_tube)
                        obj.tube(itt) = intersect(orig_tube.tube(itt), ...
                            int_polytope);
                    end
                else
                    err = SrtInvalidArgsError(['Invalid arguments, ', ...
                            'see help.']);
                        throwAsCaller(err);
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
                        throwAsCaller(SrtInvalidArgsError(['All sets in the', ...
                            'tube must be of the same dimension']));
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
        % property, otherwise the Tube object is refernced.
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
        %        https://sreachtools.github.io/license/
        % 

            if strcmp(s(1).type, '()')
                [varargout{1:nargout}] = subsref(obj.tube, s);
            else
                [varargout{1:nargout}] = builtin('subsref', obj, s); 
            end
        end

        function n_tubes = length(obj)
        % Overloading of MATLAB internal length function
        % ====================================================================
        %
        % Function to get length of tube, returns length of tube
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
        %        https://sreachtools.github.io/license/
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
        %        https://sreachtools.github.io/license/
        % 

            fprintf('Tube of %d sets\n', length(obj));
        end

        function [concat_tube_A, concat_tube_b] = concat(obj,varargin)
        % Get concatenated tube (Cartesian product of the polytopes)
        % ======================================================================
        %
        % This method computes the half-space representation of the concatenated 
        % tube. When no arguments are provided, it returns 
        % safe_set^{time_horizon} x target_set, a huge polyhedron in the
        % (obj.dim x time_horizon)-dimensional Euclidean space.
        %
        % [concat_tube_A, concat_tube_b] = concat(obj);
        %
        % The output matrices satisfy the relation that the a concatenated 
        % state vector X lies in the reach-avoid tube if and only if
        % concat_tube_A * [initial_state;X] <= concat_tube_b 
        %
        % When arguments are specified, it provides the Cartesian of
        % specific time splice mentioned. Specifically, if the half space
        % representation of sets from t=3 to 5 is desired for a given
        % tube of length 10 (sets are defined for t=0 to 9), we
        % provide
        %
        % [concat_tube_A, concat_tube_b] = concat(obj, [4 6]);
        % 
        % The +1 added is to account for MATLAB's indexing which begins
        % from 1. In other words, provide the starting and ending index of
        % interest with respect to the Tube array of polyhedrons.
        %
        % Usage: See getLowerBoundStochReachAvoid.
        %
        % ======================================================================
        %
        % [concat_tube_A, concat_tube_b] = concat(obj,varagin);
        %
        % Inputs:
        % -------
        %   time_limits - A 1x2 vector [a b] with the 1 <= a,b <= length(obj).
        %                 If a>b, then empty matrices are returned.
        %    
        % Outputs:
        % --------
        %   concat_tube_A - State matrix concatenated for tube
        %   concat_tube_b - Input matrix concatenated for tube
        %
        % Notes:
        % ------
        % * This function also serves as a delegatee for input handling.
        % 
        % ======================================================================
        % 
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://sreachtools.github.io/license/
        %
        %

            %% Construction of the concatenated tube
            tube_A_mats = cell(1, length(obj));
            [tube_A_mats{:}] = obj.tube(:).A;

            tube_b_vecs = cell(1, length(obj));
            [tube_b_vecs{:}] = obj.tube(:).b;

            %% Do we send out everything or was a slice requested?
            if nargin == 1
                % Send out everything
                concat_tube_A = blkdiag(tube_A_mats{:});
                concat_tube_b = vertcat(tube_b_vecs{:});
            else
                % Send out a slice
                time_limits = varargin{1};
                validateattributes(time_limits,{'numeric'},{'vector'});
                if length(time_limits) == 2 && min(time_limits) >= 1 &&...
                    max(time_limits) <= length(obj)
                    if time_limits(2) >= time_limits(1)
                        % a <= b => send out the appropriate slices
                        concat_tube_A =...
                            blkdiag(tube_A_mats{time_limits(1):time_limits(2)});
                        concat_tube_b =...
                            vertcat(tube_b_vecs{time_limits(1):time_limits(2)});
                    else
                        % a > b => send out empty sets
                        concat_tube_A = [];
                        concat_tube_b = [];
                    end
                else
                    % time_limits = [a,b] is not a 2x1/1x2 vector OR it does not
                    % satisfy 1<= a,b <= length(obj)
                    throw(SrtInvalidArgsError('Invalid time range'));
                end
            end
        end
        
        function [contains_flag] = contains(obj,X)
        % Check if a given concatenates state trajectory lies in the tube
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
        %      https://sreachtools.github.io/license/
        %
        %  

            if size(X, 1) ~= length(obj) * obj.dim
                throw(SrtInvalidArgsError(sprintf(['Concatenated state', ...
                    'vector/matrix provided has incorrect number of rows. ', ...
                    'Expected: %d | Got: %d'],length(obj)*obj.dim,size(X, 1))));
            end

            if size(X, 2) <= 0
                throw(SrtInvalidArgsError(['Concatenated state ', ...
                    'vector/matrix provided should have at least one column']));
            end
            
            % Matrix to store the result for each time instant, trajectory
            contains_flag_all = zeros(length(obj),size(X,2));
            % Iterate over all time indexes since MPT contains will do
            % unnecessarily n_poly x n_state comparisons
            for t_indx = 1:length(obj)
                % Parse the t^th state from X
                X_at_t_indx = X((t_indx-1)*obj.dim + 1: t_indx*obj.dim,:);
                % Use MPT's contains functionality
                contains_flag_all(t_indx,:) =...
                    obj.tube(t_indx).contains(X_at_t_indx);
            end
            % Are t^th state in the appropriate polytope for all t?
            % all for a matrix returns a row vector with checking done
            % columnwise
            contains_flag = all(contains_flag_all);
        end
        
        function v = end(obj, ~, ~)
        % Overloading of MATLAB internal end function for objects
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
        %   v - Value of end; equal to the length of the tube
        %
        % ======================================================================
        % 
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://sreachtools.github.io/license/
        %
        %
          
            v = length(obj);
        end
        
    end    
    
    methods (Static)
        function newobj = polyArray2Tube(polyArray)
        % Convert the array of polyhedra to a Tube object
        % ======================================================================
        % Inputs:
        % -------
        %   polyArray - Array of polyhedrons
        %
        % Outputs:
        % --------
        %   newobj    - Tube object
        %
        % ======================================================================
        % 
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://sreachtools.github.io/license/
        %
        %
            cellPolyArray = num2cell(polyArray);
            newobj = Tube(cellPolyArray{:});
        end
    end
end
