classdef SimpleBox
% Class definition to obtain vertices of a n-dimensional box
% ===========================================================================
%
% Class to obtain vertices of an n-dimensional box; often used for computing
% probabilities in dynamic programming recursions
%
% Usage
% -----
% % call by passing in vertices
% simpbox = SIMPLEBOX([1,1;-1,-1;-1,1;1,-1]);
% 
% % call by passing center and deltas 
% simpbox = SIMPLEBOX(0, [1, 1])
% 
% ===========================================================================
%
% SIMPLEBOX Properties:
% ---------------------
%   vertices- Array (m x n) of vertices; each vertex is a (1 x n) array
%   center  - Array (1 x n) of box center location
%   dx      - Array (1 x n) of half-lengths of box sides
%   dim     - Dimension of box (scalar)
%
% SIMPLEBOX Methods:
% ------------------
%   SimpleBox/SimpleBox        - Class constructor
%   getPolyhedron              - Get Polyhedron object for box
%   computeGaussianProbability - Compute the probability of Gaussian random 
%                                variable being in box
%
% ===========================================================================
%
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 
% 

    properties (SetAccess = private)
        vertices
        center

        % Box side half-lengths
        % =================================================================
        % 
        % Array of the half-lenghs of each side of the box
        %
        % ASCII Example of 2-d box, 'x' are vertices, 'c' is center
        %
        %  x --------------------------------- x
        %  |                 |                 |
        %  |                                   |
        %  |                dx(2)              |
        %  |                                   |
        %  |                 |                 |
        %  |                                   |
        %  | ---- dx(1) ---- c ---- dx(1) ---- |
        %  |                                   |
        %  |                 |                 |
        %  |                                   |
        %  |                dx(2)              |
        %  |                                   |
        %  |                 |                 |
        %  |                                   |
        %  x --------------------------------- x
        %
        %
        dx
        dim
    end
    
    methods
        function obj = SimpleBox(vertices, dx)
        % Class constructor for SimpleBox
        % ====================================================================
        %
        % Constructor for SimpleBox Class
        %
        % Usage:
        % ------
        % % call by passing in vertices
        % simpbox = SIMPLEBOX([1,1;-1,-1;-1,1;1,-1]);
        % 
        % % call by passing center and deltas 
        % simpbox = SIMPLEBOX(0, [1, 1])
        % 
        % ====================================================================
        %
        % obj = SimpleBox(vertices)
        % obj = SimpleBox(center, dx)
        %
        % Inputs:
        % -------
        % Vertices call form:
        %   vertices  - Array (m x n) of vertices; each vertex is a (1 x n) 
        %               array
        % 
        % Center and delta call from:
        %   center    - Array (1 x n) of box center location
        %   dx        - Array (1 x n) of half-lengths of box sides
        %
        % Outputs:
        % --------
        % obj - SimpleBox object
        %
        % ====================================================================
        %
        %   This function is part of the Stochastic Reachability Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        % 

            if size(vertices, 1) > 1
                if nargin >= 2
                    exc = SrtInvalidArgsError(['When creating a box ', ...
                        'using a vertices the delta in each dimension ', ...
                        'should not be provided; see help ', ...
                        'SimpleBox/SimpleBox']);
                    throwAsCaller(exc);
                end
                vertices = obj.sortVertices(vertices);

                [obj.center, obj.dx] = ...
                    obj.getCenterAndDeltasFromVertices(vertices);

                obj.vertices = vertices;
            else
                if nargin < 2
                    exc = SrtInvalidArgsError(['When creating a box ', ...
                        'using a center the delta in each dimension must ', ...
                        'also be provided; see help SimpleBox/SimpleBox']);
                    throwAsCaller(exc);
                end
                center = vertices;

                obj.center = center;
                obj.dx = dx;

                if length(dx) == 1
                    obj.vertices = [center - dx, center + dx]';
                elseif length(dx) == 2
                    obj.vertices = obj.get2dBoxVerticesFromCenter(center, dx);
                elseif length(dx) == 3
                    obj.vertices = obj.get3dBoxVerticesFromCenter(center, dx);
                else
                    obj.vertices = zeros(2^length(dx), size(center, 2));
                    ones_vec   = [1, -1];
                    ind_vec    = ones(1, length(dx));
                    shift_vec  = zeros(1, length(dx));
                    for i = 1:size(obj.vertices, 1)
                        for j = 1:length(dx)
                            shift_vec(j) = ones_vec(ind_vec(j));
                        end
                        obj.vertices(i, :) = center + shift_vec .* dx;

                        ind_vec(end) = ind_vec(end) + 1;
                        for j = length(dx):-1:1
                            if ind_vec(j) > 2
                                if j == 1
                                    break;
                                else
                                    ind_vec(j) = 1;
                                    ind_vec(j-1) = ind_vec(j-1) + 1;
                                end
                            else
                                break;
                            end
                        end
                    end
                end
            end
            
            obj.dim = length(obj.dx);
        end
        
        function [lb, ub] = getBounds(obj)
        % SimpleBox/getBounds  Get upper and lower bounds for simple box
        % ====================================================================
        %
        % Get upper and lower bounds for SimpleBox object.
        %
        % Usage
        % -----
        % simpbox = SimpleBox(0, [1, 1]);
        % [lb, ub] = simpbox.getBounds();
        % 
        % ====================================================================
        %
        % [lb, ub] = obj.getBounds();
        %
        % Inputs: None
        %
        % Outputs:
        % --------
        % lb - Lower bounds
        % ub - Upper bounds
        %
        % ====================================================================
        %
        %   This function is part of the Stochastic Reachability Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        % 
            lb = obj.center - obj.dx;
            ub = obj.center + obj.dx;
        end
        
        function p = computeGaussianProbability(obj, vertex_probabilites)
        % Compute the likelihood
        % for Gaussian to be in box
        % =====================================================================
        % 
        % Method to compute the likelihood for a Gaussian random variable to 
        % lie in the given box (SimpleBox object)
        %
        % Usage:
        % ------
        % simpbox = SimpleBox(0, [1, 1]);
        % p = simpbox.computeGaussianProbability(mvncdf(simpbox.vertices));
        % 
        % =====================================================================
        % 
        % p = obj.computeGaussianProbability(vertex_probabilities)
        % 
        % Inputs:
        % -------
        %   vertex_probabilities - Probabilities of Gaussian at each vertex
        %
        % Outputs:
        % --------
        %   p - Probability
        %
        % =====================================================================
        %
        %   This function is part of the Stochastic Reachability Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        % 
            
            n_dims = size(obj.vertices, 2);
            n_vertices = size(obj.vertices, 1);
            if n_vertices == 4
                p = diff(diff(reshape(vertex_probabilites, 2, [])));
            else
                p = diff(reshape(vertex_probabilites, 2, []))';
                if size(p, 1) > 1
                    box = SimpleBox(obj.center(1:n_dims-1), ...
                        obj.dx(1:n_dims-1));
                    p = box.computeGaussianProbability(p);
                end
            end
        end
        
        function poly = getPolyhedron(obj)
        % Get Polyhedron form of box
        % ================================================================
        %
        % Class method to get the MPT Polyhedron representation of the 
        % SimpleBox object
        %
        % Usage
        % -----
        % simpbox = SimpleBox(0, [1, 1]);
        % poly = simpbox.getPolyhedron();
        %
        % ================================================================
        %
        % poly = obj.getPolyhedron()
        %
        % Inputs: None
        %
        % Outputs:
        % --------
        % poly - Polyhedron object
        %
        % ================================================================
        %
        %   This function is part of the Stochastic Reachability Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        % 

            poly = Polyhedron(obj.vertices);
        end
    end
    
    methods (Static, Access = private)
        function vertices = sortVertices(vertices)
        % Sort box vertices
        % ====================================================================
        %
        % Private, static method to sort simple box vertices
        %
        % Usage:
        % ------
        % % private method of SimpleBox class
        %
        % ====================================================================
        %
        % vertices = sortVertices(vertices)
        %
        % Inputs:
        % -------
        %   vertices - m x n list of vertices, m is number of vertices and n is
        %              dimension of system
        %
        % Outputs:
        % --------
        %   vertices - m x n list of vertices, m is number of vertices and n is
        %              dimension of system
        %
        % ====================================================================
        %
        %   This function is part of the Stochastic Reachability Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        % 
        
            for i = size(vertices, 2):-1:1
                vertices = sortrows(vertices, i);
            end
        end

        function [center, dx] = getCenterAndDeltasFromVertices(vertices)
            dx = zeros(1, size(vertices, 2));

            n_verts = size(vertices, 1);
            for i = 1:size(vertices, 2)
                vert_diff = diff(sortrows(vertices, i, 'ascend'));
                dx(i) = vert_diff(n_verts/2, i) / 2;
            end

            vertices = SimpleBox.sortVertices(vertices);
            center = vertices(1, :) + dx;
        end

        function vertices = get2dBoxVerticesFromCenter(center, dx)
        % SReachTools/SimpleBox/get2dBoxVerticesFromCenter
        % ====================================================================
        %
        % Priavte, static method to compute points that create a bounding box 
        % around a given ceneter and spacing in each direction for 2-dimensional 
        % systems. The small system can be hard-coded for increased speed.
        % 
        % Usage:
        % ------
        % % privae method of SimpleBox
        %
        % ====================================================================
        %
        % Inputs:
        % -------
        %   center - Current SpaceGrid grid point
        %   dx     - State SpaceGrid deltas
        %
        % Outputs:
        % --------
        %   vertices - Array of vertices
        %
        % ====================================================================
        %
        %   This function is part of the Stochastic Reachability Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        % 

            vertices = zeros(4, length(dx));
            vertices(1, :) = center + [1, 1] .* dx;
            vertices(2, :) = center + [1, -1] .* dx;
            vertices(3, :) = center + [-1, 1] .* dx;
            vertices(4, :) = center + [-1, -1] .* dx;
            
            vertices = SimpleBox.sortVertices(vertices);
        end

        function vertices = get3dBoxVerticesFromCenter(center, dx)
        % SReachTools/SimpleBox/get3dBoxVerticesFromCenter
        % ====================================================================
        %
        % Private, static method to compute points that create a bounding box 
        % around a given ceneter and spacing in each direction for 2-dimal 
        % systems. The small system can be hard-coded for increased speed.
        % 
        % Usage:
        % ------
        % % private method of SimpleBox
        %
        % ====================================================================
        %
        % Inputs:
        % -------
        %   center - Current SpaceGrid grid point
        %   dx     - State SpaceGrid deltas
        %
        % Outputs:
        % --------
        %   vertices - Array of vertices
        % 
        % ====================================================================
        %
        %   This function is part of the Stochastic Reachability Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        % 

            vertices = zeros(8, length(dx));
            vertices(1, :) = center + [1,   1,  1] .* dx;
            vertices(2, :) = center + [1,   1, -1] .* dx;
            vertices(3, :) = center + [1,  -1,  1] .* dx;
            vertices(4, :) = center + [1,  -1, -1] .* dx;
            vertices(5, :) = center + [-1,  1,  1] .* dx;
            vertices(6, :) = center + [-1,  1, -1] .* dx;
            vertices(7, :) = center + [-1, -1,  1] .* dx;
            vertices(8, :) = center + [-1, -1, -1] .* dx;
            
            vertices = SimpleBox.sortVertices(vertices);
        end
    end
end
