classdef SpaceGrid
% SReachTools/SpaceGrid  Create a state space grid object
% =============================================================================
%
% Class to hold the gridding of a particular space, e.g. state or input.
%
% Usage:
% ------
% % Define a 2-dimensional state-space grid that extends from x = [-1, 1],
% % y = [-1, 1] with 100 points in each dimension
%
% grid = SpaceGrid([-1, -1], [1, 1], 100)
%
% % Can also define different dimensional spacings
%
% grid = SpaceGrid([-1, -1], [1, 1], [100, 50])
% 
% =============================================================================
%
% SpaceGrid Properties:
% ---------------------
%   grid         - Array of grid vectors, size prod(n_points) x dim
%   lower_bounds - Lower bounds provided during construction
%   upper_bounds - Upper bounds provided during construction
%   n_points     - Number of points in grid in each dimension
%   grid_delta   - Grid spacing, spacing between two grid points is 
%                  2*grid_delta(i); 'i' being the dimension of interest
%   dim          - Total number of dimensions in the grid
% 
% SpaceGrid Methods:
% ------------------
%   SpaceGrid/SpaceGrid      - Constructor
%   getIndicatorVectorForSet - Method to get an indicator vector of which grid 
%                              points are in a Polyhedron set
%   getMeshGrids             - Get the associated MATLAB mesh grids for the 
%                              space; only works for 2 or 3-dimensional systems
%   getExternalGrid          - Method to get an external grid for the current
%                              grid
%   plotGridProbability      - Helper method for plotting dynamic programming
%                              probabilities on the grid
% 
% =============================================================================
%
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 
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
        actual_lb
        actual_ub
        total_points
        is_external
    end
    methods
        function obj = SpaceGrid(lb, ub, n_points, external_flag)
        % SReachTools/SpaceGrid/SpaceGrid  Constructor
        % ====================================================================
        %
        % SpaceGrid class constructor
        %
        % Usage:
        % ------
        % % Define a 2-dimensional state-space grid that extends from 
        % % x = [-1, 1], y = [-1, 1] with 100 points in each dimension
        %
        % grid = SpaceGrid([-1, -1], [1, 1], 100)
        %
        % % Can also define different dimensional spacings
        %
        % grid = SpaceGrid([-1, -1], [1, 1], [100, 50])
        %   
        % ====================================================================
        %
        % Inputs:
        % -------
        %   lb            - Lower bounds
        %   ub            - Upper bounds
        %   n_points      - Number of points in each grid dimension
        %   external_flag - (optional) External Flag argument
        %
        % Outputs:
        % --------
        %   obj - SpaceGrid object
        %
        % ====================================================================
        %
        %   This function is part of the Stochastic Reachability Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        % 

            if nargin < 4
                external_flag = 'internal';
                is_external = true;
            end
            
            validateattributes(external_flag, {'char', 'string'}, {'nonempty'});
            
            external_flag = lower(external_flag);
            
            if strcmp(external_flag, 'internal')
                obj.is_external = false;
            elseif strcmp(external_flag, 'external')
                obj.is_external = true;
            else
                exc = SrtInvalidArgsError(['Flag for external ', ...
                        'grid creation must be either ''internal'' or ', ...
                        '''external''']);
                throwAsCaller(exc);
            end
            
            % validate that the lower and upper bounds are positie integer
            % values
            validateattributes(lb, {'numeric'}, {'vector'})
            validateattributes(ub, {'numeric'}, {'vector'})
            
            % check if lower and upper are same dimension
            if length(lb) ~= length(ub)
                exc = SrtInvalidArgsError(['Lower and upper bounds ', ...
                    'must be equivalent in length (dimension).']);
                throwAsCaller(exc);
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
                    exc = SrtInvalidArgsError(['Number of points ', ...
                        'must be either a scalar or a vector of ', ...
                        'equivalent length of the lower and upper bounds.']);
                    throwAsCaller(exc);
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
                obj.grid_delta(i) = (ub(i) - lb(i)) / (2 * (n_points(i) - 1));
            end
            
            % depending on if the grid is internal or external need to set the
            % actual lower and upper bounds
            switch(external_flag)
                case 'internal'
                    obj.actual_lb = lb;
                    obj.actual_ub = ub;
                case 'external'
                    obj.actual_lb = lb + obj.grid_delta;
                    obj.actual_ub = ub - obj.grid_delta;
                otherwise
                    exc = SrtInvalidArgsError(['Flag for external ', ...
                        'grid creation must be either ''internal'' or ', ...
                        '''external'''])
                    throwAsCaller(exc);
            end
            
            
            % start making the grid
            % initialize zeros
            obj.grid = zeros(obj.total_points, length(lb));
            inds = ones(1, length(n_points));
            for i = 1:obj.total_points
                obj.grid(i, :) = obj.getGridVectorFromInds(inds);
                
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
        % SReachTools/SpaceGrid/getIndicatorVectorForSet  Get indicator vector
        % ====================================================================
        % 
        % Get indicator vector for the grid points which lie in a Polyhedron
        % set, s.
        %
        % Usage:
        % ------
        %   grid = SpaceGrid([-1, -1], [1, 1], 100);
        %   ind_vector = grid.getIndicatorVectorForSet(Polyhedron(...
        %       'lb', [0,0], 'ub', [1, 1]));
        %
        % ====================================================================
        % 
        % ind_vector = getIndicatorVectorForSet(obj, s)
        %
        % Inputs:
        % -------
        %   s - Polyhedron object set
        % 
        % Outputs:
        % --------
        %   ind_vector - Indicator vector for grid points in set s
        %
        % ====================================================================
        %
        %   This function is part of the Stochastic Reachability Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        % 

            validateattributes(s, {'Polyhedron'}, {'nonempty'})
            ind_vector = s.contains(obj.grid');
        end

        function sortGrid(obj)
        % SReachTools/SpaceGrid/sortGrid  Sort space grid vectors
        % ====================================================================
        % 
        % Sort space grid vectors, ascending
        %
        % Usage:
        % ------
        %   grid = SpaceGrid([-1, -1], [1, 1], 100);
        %   grid.sortGrid();
        %
        % ====================================================================
        %
        % sortGrid(obj)
        %
        % Inputs:  None
        % Outputs: None
        % 
        % ====================================================================
        %
        %   This function is part of the Stochastic Reachability Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        %
        
            for i = size(obj.grid, 2):-1:1
                obj.grid = sortrows(obj.grid, i);
            end
        end
        
        function varargout = getMeshGrids(obj)
        % SReachTools/SpaceGrid/getMeshGrids  Get MATLAB meshgrids
        % ====================================================================
        % 
        % Get MATLAB meshgrids for the SpaceGrid object; only works for grids 
        % that are 2 or 3-dimensional
        % 
        % Usage:
        % ------
        %   grid = SpaceGrid([-1, -1], [1, 1], 100);
        %   [X,Y] = grid.getMeshGrids();
        % 
        % ====================================================================
        % 
        % [X,Y] = obj.getMeshGrids();
        % [X,Y,Z] = obj.getMeshGrids();
        %
        % Inputs: None
        % 
        % Outputs:
        % --------
        %   X - x-meshgrid
        %   Y - y-meshgrid
        %   Z - z-meshgrid
        %
        % ====================================================================
        %
        %   This function is part of the Stochastic Reachability Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        %

            if obj.dim > 3
                exc = SrtInternalError(['Can only get meshgrid for 2 ', ...
                    'or 3-dimentional systems']);
                throw(exc);
            end
            
            if obj.dim == 2
                [mesh_x, mesh_y] = meshgrid(...
                    obj.actual_lb(1):2*obj.grid_delta(1):obj.actual_ub(1), ...
                    obj.actual_lb(2):2*obj.grid_delta(2):obj.actual_ub(2));
                    
                varargout = {mesh_x, mesh_y};
            else
                [mesh_x, mesh_y, mesh_z] = meshgrid(...
                    obj.actual_lb(1):2*obj.grid_delta(1):obj.actual_ub(1), ...
                    obj.actual_lb(2):2*obj.grid_delta(2):obj.actual_ub(2), ...
                    obj.actual_lb(3):2*obj.grid_delta(3):obj.actual_ub(3));
                    
                varargout = {mesh_x, mesh_y, mesh_z};
            end
        end
        
        function ext_grid = getExternalGrid(obj)
            if ~obj.is_external
                ext_grid = SpaceGrid(...
                    obj.lower_bounds, ...
                    obj.upper_bounds, ...
                    obj.n_points + 1, ...
                    'external');
            else
                exc = SrtInternalError(['Cannot create an external ', ...
                    'grid from an external grid']);
                throw(exc);
            end
        end
        
        function plotGridProbability(obj, grid_probability)
        % SReachTools/SpaceGrid/plotGridProbability  Plot grid probability
        % ====================================================================
        % 
        % Perform surface plot of 2-dimensional grid probability
        % 
        % Usage:
        %   grid = SpaceGrid([-1, -1], [1, 1], 100);
        %   grid.plotGridProbability(mvncdf(grid.grid));
        % 
        % ====================================================================
        % 
        % obj.plotGridProbability
        % 
        % Inputs:  None
        % Outputs: None
        % 
        % ====================================================================
        %
        %   This function is part of the Stochastic Reachability Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        %
            if obj.dim > 2
                exc = SrtInternalError(['Can only plot for 1 or ', ...
                    '2-dimentional systems']);
                throw(exc);
            end
            
            if obj.dim == 2
                [X,Y] = obj.getMeshGrids();
                surf(X, Y, reshape(grid_probability, obj.n_points));
            else
                % something here later...    
            end
            
        end
    end
    
    methods (Hidden)
        function plotGrid(obj)
        % SReachTools/SpaceGrid/plotGridProbability  Plot grid 
        % ====================================================================
        % 
        % Hidden method to plot 2 or 3-dimensional grid
        % 
        % Usage:
        %   grid = SpaceGrid([-1, -1], [1, 1], 100);
        %   grid.plotGrid();
        % 
        % ====================================================================
        % 
        % obj.plotGrid
        % 
        % Inputs:  None
        % Outputs: None
        % 
        % ====================================================================
        %
        %   This function is part of the Stochastic Reachability Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        %

            if obj.dim > 3
                exc = SrtInternalError(['Can only plot the grid for ', ...
                    '2 or 3-dimensional grids'])
                throw(exc);
            end
            
            if obj.dim == 2
                [X,Y] = obj.getMeshGrids();
                X = reshape(X, [], 1);
                Y = reshape(Y, [], 1);
                
                scatter(X, Y);
            else
                [X,Y,Z] = obj.getMeshGrids();
                X = reshape(X, [], 1);
                Y = reshape(Y, [], 1);
                Z = reshape(Z, [], 1);
                
                scatter3(X, Y, Z);
            end
        end
        
        function arrays = getGridArrays(obj)
        % SReachTools/SpaceGrid/getGridArrays  Get arrays for grid vectors for each
        % dimension
        % ====================================================================
        % 
        % Hidden method to get the grid vectors that make up the grid points
        % independently for each dimension
        % 
        % Usage:
        % ------
        %   grid = SpaceGrid([-1, -1], [1, 1], 100);
        %   vecs = grid.getGridsArarys();
        % 
        % ====================================================================
        % 
        % arrays = obj.getGridArrays()
        %
        % Inputs: None
        % Outputs:
        % --------
        %   arrays - Cell array of grid arrays
        % 
        % ====================================================================
        %
        %   This function is part of the Stochastic Reachability Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        %

            arrays = cell(1, obj.dim);
            
            for i = 1:obj.dim
                arrays{i} = obj.lower_bounds(i)+obj.grid_delta(i):...
                    2*obj.grid_delta(i):...
                    obj.upper_bounds(i)-obj.grid_delta(i);
            end
        end
    end
    
    methods (Access = private)
        function grid_vector = getGridVectorFromInds(obj, inds)
        % SReachTools/SpaceGrid/getGridVectorFromInds  Get the grid vector from 
        % indices
        % ====================================================================
        %
        % Private method, get the grid vector from indices for each value in
        % its indipendent grid array
        % 
        % Usage: Private method, thus not callable from command prompt
        %
        % ====================================================================
        %
        % obj.getGridVectorFromInds(inds)
        % 
        % Inputs:
        % -------
        %   inds - Indices for each point in the individual grid array
        % 
        % Outputs:
        % --------
        %   grid_vector - Grid vector
        %
        % ====================================================================
        %
        %   This function is part of the Stochastic Reachability Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        %
        
            % validate inputs
            validateattributes(inds, {'numeric'}, ...
                {'>', 0, 'integer', 'vector'});
            
            % make sure that the length of the indices matches the bounds
            % dimension
            if length(inds) ~= length(obj.lower_bounds)
                exc = SrtInternalError(['Length/dimension of indices ', ...
                    'do not match the length/dimensions of the bounds']);
                throw(exc);
            end
            
            % get the grid vector
            grid_vector = zeros(1, length(inds));
            for i = 1:length(inds)
                % line calculation
                a = obj.actual_ub(i);
                b = obj.actual_lb(i);
                N = obj.n_points(i);
                grid_vector(i) = ((a - b) / (N - 1)) * (inds(i) - 1) + b;
            end
        end
    end
        
end