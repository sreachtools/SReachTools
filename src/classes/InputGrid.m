classdef InputGrid
% SReachTools/InputGrid: Create a input space grid object
% ============================================================================
%
% Defines an input space grid used for dynamic programming computations
%
% Usage:
% ------
%
% % Define a 2-dimensional state-space grid that extends from x = [-1, 1],
% % y = [-1, 1] with 100 points in each dimension
%
% grid = InputGrid([-1, -1], [1, 1], 100)
%
% % Can also define different dimensional spacings
%
% grid = InputGrid([-1, -1], [1, 1], [100, 50])
%   
% ============================================================================
%
% InputGrid Properties:
% ---------------------
%   grid         - Array of grid vectors, size prod(n_points) x dim
%   lower_bounds - Lower bounds provided during construction
%   upper_bounds - Upper bounds provided during construction
%   n_points     - Number of points in grid in each dimension
%   grid_delta   - Grid spacing, spacing between two grid points is 
%                  2*grid_delta(i); 'i' being the dimension of interest
%   dim          - Total number of dimensions in the grid
% 
% InputGrid Methods:
% ------------------
%   InputGrid/InputGrid      - Class constructor
%   getIndicatorVectorForSet - Method to get an indicator vector of which grid 
%                              points are in a Polyhedron set
%
% Notes:
%   - After the 'external' option is removed from the SpaceGrid class there will
%     be very little difference between InputGrid and SpaceGrid. The classes
%     will be merged for convenience.
% 
% ============================================================================
%
%   This function is part of the Stochastic Optimal Control Toolbox.
%   License for the use of this function is given in
%        https://github.com/abyvinod/SReachTools/blob/master/LICENSE
% 

    properties (SetAccess = private)
        grid
        lower_bounds
        upper_bounds
        n_points
        grid_delta
        dim
    end
    properties (Access = private)
        total_points
    end
    methods
        function obj = InputGrid(lb, ub, n_points)
        % SReachTools/InputGrid/InputGrid
        % ====================================================================
        %
        % InputGrid class constructor
        %
        % Usage:
        % ------
        % % Define a 2-dimensional state-space grid that extends from x = 
        % % [-1, 1], y = [-1, 1] with 100 points in each dimension
        %
        % grid = InputGrid([-1, -1], [1, 1], 100)
        %
        % % Can also define different dimensional spacings
        %
        % grid = InputGrid([-1, -1], [1, 1], [100, 50])
        % 
        % ====================================================================
        %
        % obj = InputGrid(lb, up, n_points)        
        % 
        % Inputs:
        % -------
        %   lb       - Lower bounds in each dimension
        %   ub       - Upper bounds in each dimension
        %   n_points - Number of point for each dimension
        % 
        % Outputs:
        % --------
        %   obj - InputGrid object
        % 
        % ====================================================================
        %
        %   This function is part of the Stochastic Optimal Control Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/abyvinod/SReachTools/blob/master/LICENSE
        %   

            % validate that the lower and upper bounds are positie integer
            % values
            validateattributes(lb, {'numeric'}, {'vector'})
            validateattributes(ub, {'numeric'}, {'vector'})
            
            % check if lower and upper are same dimension
            if length(lb) ~= length(ub)
                error('SReachTools:invalidArgs', ['Lower and upper bounds must ', ...
                    'be equivalent in length (dimension).']);
            end
            
            % if number of points is a scalar then copy amount for each
            % dimension
            validateattributes(n_points, {'numeric'}, ...
                {'integer', '>', 0, 'vector'})
            if length(n_points) == 1
                n_points = n_points * ones(size(lb));
            else
                % need to check to ensure that the numer of points is the same
                % as the number of dimensions of the lower and upper bounds
                if length(n_points) ~= length(lb)
                    error('SReachTools:invalidArgs', ['Number of points must ', ...
                        'be either a scalar or a vector of equivalent ', ...
                        'length of the lower and upper bounds.']);
                end
            end
            
            % set object properties
            obj.lower_bounds = lb;
            obj.upper_bounds = ub;
            obj.dim          = length(lb);
            obj.n_points     = n_points;
            obj.total_points = prod(n_points);
            
            % compute the grid spacing delta
            obj.grid_delta = zeros(1, length(lb));
            for i = 1:length(obj.grid_delta)
                obj.grid_delta(i) = (ub(i) - lb(i)) / (n_points(i) - 1);
            end
            
            % start making the grid
            % initialize zeros
            obj.grid = zeros(obj.total_points, length(lb));
            inds = ones(1, length(n_points));
            for i = 1:obj.total_points
                obj.grid(i,:) = obj.getGridVectorFromInds(inds);

                inds(end) = inds(end) + 1;
                for j = length(inds):-1:1
                    if inds(j) > obj.n_points(j)
                        if j == 1
                            break;
                        else
                            inds(j) = 1;
                            inds(j-1) = inds(j-1) + 1;
                        end
                    else
                        break;
                    end
                end
            end
        end
        
        function ind_vector = getIndicatorVectorForSet(obj, s)
        % SReachTools/InputGrid/getIndicatorVectorForSet Get indicator vector for
        % points in grid that lie in set 
        % ====================================================================
        %
        % Class method returning an indicator vector (vector of zeros and ones)
        % for all points in the grid that lie in a Polyhedral set
        %
        % Usage:
        % ------
        % grid = InputGrid([-1, -1], [1, 1], 100);
        % s = Polyhedron('lb', [0, 0], 'ub', [1, 1]);
        % vec = grid.getIndicatorVectorForSet(s);
        % 
        % ====================================================================
        %
        % ind_vector = getIndicatorVectorForSet(obj, s)       
        % 
        % Inputs:
        % -------
        %   obj - InputGrid object
        %   s   - Polyhedron set
        % 
        % Outputs:
        % --------
        %   ind_vector - Indicator vector (n_points x 1)
        % 
        % ====================================================================
        %
        %   This function is part of the Stochastic Optimal Control Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/abyvinod/SReachTools/blob/master/LICENSE
        % 

            validateattributes(s, {'Polyhedron'}, {'nonempty'})
            ind_vector = s.contains(obj.grid');
        end
    end
    
    methods (Hidden)
        function grid_vector = getGridVectorFromInds(obj, inds)
        % SReachTools/InputGrid/getGridVectorFromInds  Get grid vector
        % ====================================================================
        %
        % Class method returning an indicator vector (vector of zeros and ones)
        % for all points in the grid that lie in a Polyhedral set
        %
        % Usage: Hidden method
        % 
        % ====================================================================
        %
        % grid_vector = getGridVectorFromInds(obj, inds)    
        % 
        % Inputs:
        % -------
        %   obj  - InputGrid object
        %   inds - ???
        % 
        % Outputs:
        % --------
        %   grid_vector - ???
        % 
        % ====================================================================
        %
        %   This function is part of the Stochastic Optimal Control Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/abyvinod/SReachTools/blob/master/LICENSE
        %    

            % validate inputs
            validateattributes(inds, {'numeric'}, ...
                {'>', 0, 'integer', 'vector'});
            
            % make sure that the length of the indices matches the bounds
            % dimension
            if length(inds) ~= length(obj.lower_bounds)
                error('SReachTools:internal', ['Length/dimension of indices do ', ...
                    'not match the length/dimensions of the bounds']);
            end
            
            % get the grid vector
            grid_vector = zeros(1, length(inds));
            for i = 1:length(inds)
                % line calculation
                a = obj.upper_bounds(i);
                b = obj.lower_bounds(i);
                N = obj.n_points(i);
                grid_vector(i) = ((a - b) / (N - 1)) * (inds(i) - 1) + b;
            end
        end
    end
end