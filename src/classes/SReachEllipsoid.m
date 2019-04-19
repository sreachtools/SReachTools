classdef SReachEllipsoid
% Creates an ellipsoid (x - c)^T Q^{-1} (x-c) <= 1
% =============================================================================
% 
% Create an ellipsoid object defined by the equation
% 
%   E = { x \in R^{n} : (x - c)^{T} Q^{-1} (x - c) <= 1 }
% 
% These ellipsoid objects are often used when creating bounded disturbances
% for the SReachSet Lagrangian methods when using Gaussian disturbances
% 
% Usage:
% ------
% % Create unit ellipsoid
% sre = SReachEllipsoid([0;0], eye(2));
% 
% =============================================================================
%
% SReachEllipsoid Properties:
% ---------------------------
%   center                          - Center of the ellipsoid (c)
%   shape_matrix                    - Shape matrix of the ellipsoid Q
%   dim                             - Dimension of the ellipsoid
% 
% SReachEllipsoid Methods:
% ------------------------
%   SReachEllipsoid/SReachEllipsoid - Constructor
%   support                         - Support function of the ellipsoid
%   contains                        - Checks if a point (column vector) or a
%                                     collection of points (matrix of column
%                                     vectors) is within the ellipsoid
%
% Apart from these methods, the following commands work
%   disp                            - Displays critical info about the ellipsoid
%   F * ell | ell * F               - Multiplication of ellipsoid by a n x dim -
%                                     dimensional matrix or a scalar F
%   F + ell | ell + F | ell + poly  - Add a deterministic vector/scalar to an
%                                     ellipsoid | Overapproximate the Minkowski
%                                     sum of a polyhedron with ellipsoid
% 
% Notes:
% ------
% * The ellipsoid can be full-dimensional (Q non-singular) or be a 
%   lower-dimensional ellipsoid embedded in a high dimensional space (Q 
%   singular)
%
% ===========================================================================
%
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://sreachtools.github.io/license/
% 
% 

    properties
        % SReachEllipsoid/center
        % ==================================================================
        % 
        % Column vector indicating the center of the ellipsoid
        %
        center

        % SReachEllipsoid/shape_matrix
        % ==================================================================
        % 
        % Shape matrix for the ellipsoid
        % 
        shape_matrix

        % SReachEllipsoid/dim
        % ==================================================================
        % 
        % Dimension of the ellipsoid dimension
        % 
        dim
    end
    methods
        function obj = SReachEllipsoid(center, shape_matrix)
        %  Constructor for SReachEllipsoid class
        % ====================================================================
        %
        % Inputs:
        % -------
        %   center       - Center of the ellipsoid; must be a column vector
        %   shape_matrix - Shape matrix of the ellipsoid
        %
        % Outputs:
        % --------
        %   obj          - SReachEllipsoid object
        %
        % =====================================================================
        % 
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://sreachtools.github.io/license/
        % 
        % 
        
            % Input parsing
            inpar = inputParser();
            inpar.addRequired('center', @(x) validateattributes(x,...
                {'numeric'}, {'column','nonempty'}));
            inpar.addRequired('shape_matrix', @(x) validateattributes(x,...
                {'numeric'}, {'square','nonempty'}));

            try
                inpar.parse(center, shape_matrix);
            catch err
                exc = SrtInvalidArgsError.withFunctionName();
                exc = exc.addCause(err);
                throwAsCaller(exc);
            end
            
            obj.center = center;
            obj.shape_matrix = shape_matrix;
            obj.dim = size(obj.center,1);

            % Check if the center and shape_matrix are of correct dimensions
            if size(obj.center, 1) ~= size(obj.shape_matrix, 1)
                throwAsCaller(SrtInvalidArgsError(['Center and shape matrix',...
                    ' have different dimensions']));
            end
            if ~issymmetric(shape_matrix)
                % Compute the symmetric component of it
                symm_shape_matrix = (shape_matrix+shape_matrix')/2;
                % Max error element-wise
                max_err = max(max(abs(shape_matrix - symm_shape_matrix)));
                if max_err > eps
                    warning('SReachTools:runtime',sprintf(['Non-symmetric ', ...
                        'shape matrix made symmetric (max element-wise',...
                        ' error: %1.3e)!'], max_err));
                end
                obj.shape_matrix = symm_shape_matrix;
            end
            % For some reason, -eps alone is not enough?
            min_eig_val = min(eig(obj.shape_matrix));
            if  min_eig_val <= -2*eps
                throwAsCaller(SrtInvalidArgsError(['Covariance ',...
                    'matrix can not have negative eigenvalues']));
            elseif min_eig_val <= eps
                warning('SReachTools:runtime',['Creating an',...
                    ' SReachEllipsoid, which might be lower-dimensional']);
            end
        end

        
        function disp(obj)
        % Override of MATLAB internal display
        % ====================================================================
        % 
        % Overriding of MATLAB built-in display function for the class
        %
        % ====================================================================
        % 
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://sreachtools.github.io/license/
        % 
        %
            
            fprintf('%d-dimensional ellipsoid\n', obj.dim);
        end
        
        function val = support(obj, l)
        % Support function of the ellipsoid object
        % ====================================================================
        %
        % Inputs:
        % -------
        %   l   - A query column vector or a collection of query vectors stacked 
        %         as columns
        %
        % Outputs:
        % --------
        %   val - max_{y \in ellipsoid} l'*y
        %
        % =====================================================================
        % 
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://sreachtools.github.io/license/
        % 
        % 
        
            if size(l,1) ~= obj.dim
                throwAsCaller(SrtInvalidArgsError('l has incorrect dimensions.'));
            end
            % cholesky > obj.shape_matrix = sqrt_shape_matrix'*sqrt_shape_matrix
            [sqrt_shape_matrix, p] = chol(obj.shape_matrix);    
            if p > 0
                % Non-positive definite matrix can not use Cholesky's decompose
                % Use sqrt to obtain a symmetric non-sparse square-root matrix
                sqrt_shape_matrix = sqrt(obj.shape_matrix);
            end
            % Hence, we need the transpose
            val = l'* obj.center + norms(l'*sqrt_shape_matrix', 2, 2);
        end
        
        function newobj=mtimes(obj, F)
        % Override of MATLAB multiplication command
        % ====================================================================
        % 
        % Inputs:
        % -------
        %   obj - SReachEllipsoid object
        %   F   - Linear transformation matrix for multiplication
        %
        % Outputs:
        % --------
        %   newobj - SReachEllipsoid object (F*obj)
        %
        % ====================================================================
        % 
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://sreachtools.github.io/license/
        % 
        %
            
            switch [class(obj), class(F)]
                case ['SReachEllipsoid','double']
                    % All ok
                case ['double', 'SReachEllipsoid']
                    % Need to switch the arguments
                    Ftemp = obj;
                    obj = F;
                    F = Ftemp;
                otherwise
                    throwAsCaller(SrtInvalidArgsError(sprintf(['Operation *',...
                       ' not defined between *%s, %s'], class(obj), class(F))));
            end
            newobj=SReachEllipsoid(F * obj.center, F * obj.shape_matrix * F');            
        end
        
        function newobj = plus(obj, v)
        % Override of MATLAB plus command
        % ====================================================================
        % 
        % Inputs:
        % -------
        %   obj - Ellipsoid object
        %   v   - Deterministic vector to be added to the random vector OR
        %         a Polytope object
        %
        % Outputs:
        % --------
        %   newobj - Ellipsoid obj (obj + v) for deterministic vector/scalar v
        %            Polyhedron obj (obj \oplus v) for polytopic v (overapprox)
        %
        % Notes:
        % ------
        % * For a polytopic v, newobj is an (Polyhedron overapproximation of the
        %   minkowski sum, computed via sampling the support function.
        %
        % ====================================================================
        % 
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://sreachtools.github.io/license/
        % 
        %
            
            summands_type = [];
            switch [class(obj), class(v)]
                case ['SReachEllipsoid','double']
                    % Check dimensions
                    if ~isequal(size(v), [obj.dim 1]) && ~isequal(size(v),[1 1])
                        throwAsCaller(SrtInvalidArgsError(['Mismatch in ',...
                            'dimensions of the SReachEllipsoid and v']));
                    end
                    % Set the flag for the type of summation
                    summands_type = 'determ_vec_plus_ell';
                case ['double', 'SReachEllipsoid']
                    % Need to switch the arguments
                    vtemp = obj;
                    obj = v;
                    v = vtemp;
                    % Check dimensions
                    if ~isequal(size(v), [obj.dim 1]) && ~isequal(size(v),[1 1]) 
                        throwAsCaller(SrtInvalidArgsError(['Mismatch in ',...
                            'dimensions of the random vector and v']));
                    end
                    % Set the flag for the type of summation
                    summands_type = 'determ_vec_plus_ell';
                case ['SReachEllipsoid','Polyhedron']                
                    % Check dimensions
                    if v.Dim ~= obj.dim
                        throwAsCaller(SrtInvalidArgsError(['Mismatch in ',...
                            'dimensions of the SReachEllipsoid and ',...
                            'Polyhedron v']));
                    end
                    % Set the flag for the type of summation
                    summands_type = 'ell_plus_poly';
                case ['Polyhedron','SReachEllipsoid']                
                    % Need to switch the arguments
                    vtemp = obj;
                    obj = v;
                    v = vtemp;
                    % Check dimensions
                    if v.Dim ~= obj.dim
                        throwAsCaller(SrtInvalidArgsError(['Mismatch in ',...
                            'dimensions of the SReachEllipsoid and ',...
                            'Polyhedron v']));
                    end
                    % Set the flag for the type of summation
                    summands_type = 'ell_plus_poly';
                otherwise
                    throwAsCaller(SrtInvalidArgsError(sprintf(['Operation +',...
                       ' not defined between %s and %s'], class(obj),...
                       class(v))));
            end
            switch summands_type
                case 'determ_vec_plus_ell'
                    if isequal(size(v),[1 1])
                        % F is a scalar
                        v = v * ones(obj.dim,1);
                    end
                    newobj = SReachEllipsoid(obj.center + v, obj.shape_matrix);
                case 'ell_plus_poly'
                    new_b = v.support(v.A') + obj.support(v.A');
                    newobj = Polyhedron('H',[v.A new_b]);
                otherwise
                    % Will never come here
            end
        end
        
        function flag = contains(obj, test_points)
        % Checks if a point (column vector) or a collection of points (matrix of
        % column vectors) is within the ellipsoid
        % ====================================================================
        % 
        % Inputs:
        % -------
        %   obj         - Ellipsoid object
        %   test_points - Point (column vector) or a N collection of points
        %                 (matrix of column vectors) is within the ellipsoid
        %
        % Outputs:
        % --------
        %   newobj      - Boolean vector Nx1 that describe the containment
        %
        % Notes:
        % ------
        % * Requires CVX for vectorized norm.
        %
        % ====================================================================
        % 
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://sreachtools.github.io/license/
        % 
        %
            if size(test_points, 1) ~= obj.dim
                throwAsCaller(SrtInvalidArgsError(['test_points must be a ',...
                    'matrix with obj.dim-dimensional column vectors']));
            end
            centered_test_points = test_points - repmat(obj.center, 1,...
                size(test_points, 2));
            % Non-positive definite matrix can not use Cholesky's decompose
            % Use sqrt to obtain a symmetric non-sparse square-root matrix
            inv_shape_matrix_sqrt = sqrt(inv(obj.shape_matrix));
            % Hence, we need the transpose
            flag = norms(centered_test_points' * inv_shape_matrix_sqrt',2,2)<=1;
        end
    end
end
