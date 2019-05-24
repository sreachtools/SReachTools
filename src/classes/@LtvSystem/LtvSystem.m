classdef LtvSystem
% Create a discrete-time LTV system object
% ============================================================================
%
% Defines a discrete-time LTV system that is:
%     - control-free and disturbance-free, or
%     - controlled but disturbance-free, or
%     - perturbed (stochastic/uncertain) but control-free, or
%     - controlled and perturbed (stochastic/uncertain).
%
% Perturbation can be either:
%     - a bounded uncertainity with no stochastic information
%     - a RandomVector object
%
%  Usage:
%  ------
%  % Define a double integrator system:
%
%  T = 0.5;
%  sys = LtvSystem('StateMatrix', [1, T; 0, 1], ...
%                  'InputMatrix', [T^2/2;T], ...
%                  'InputSpace', Polyhedron('lb', -1, 'ub', 1), ...
%                  'DisturbanceMatrix', [T^2/2;T], ...
%                  'Disturbance', Polyhedron('lb', -1, 'ub', 1));
%   
% ============================================================================
%
% LtvSystem Properties:
% ---------------------
%   state_mat       - State matrix (Square matrix, state_dim x state_dim)
%   input_mat       - Input matrix (Matrix, state_dim x input_dim)
%   input_space     - Input space (empty / Polyhedron)
%   dist_mat        - Disturbance matrix (Matrix, state_dim x dist_dim)
%   dist            - Disturbance object (empty/Polyhedron/RandomVector)     
%   state_dim       - State dimension (scalar)   
%   input_dim       - Input dimension (scalar)  
%   dist_dim        - Disturbance dimension (scalar)
% 
% LtvSystem Methods:
% ------------------
%   LtvSystem/LtvSystem   - Constructor
%   getConcatInputSpace   - Get concatenated input space
%   getConcatMats         - Get concatenated state, input, and disturbance
%                           matrices
%   islti                 - Get logical value 1 if system is LTI
%   isltv                 - Get logical value 1 if system is LTV (strictly)
% 
% Notes:
% ------
% * EXTERNAL DEPENDENCY: Uses MPT3 to define input,robust disturbance space
%
% =============================================================================
%
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://sreachtools.github.io/license/
%
%

    properties (SetAccess = protected)
        % LTVSystem/state_mat
        % ====================================================================
        % 
        % System state matrix
        %   - MATLAB matrix for LTI system
        %   - Anonymous functions for LTV system, maps time index to matrix
        % 
        state_mat

        % LTVSystem/state_dim
        % ====================================================================
        % 
        % Dimension of the state vector
        % 
        state_dim

        % LTVSystem/input_mat
        % ====================================================================
        % 
        % System input matrix
        %   - MATLAB matrix for LTI system
        %   - Anonymous functions for LTV system, maps time index to matrix
        % 
        input_mat

        % LTVSystem/input_space
        % ====================================================================
        % 
        % System input space, polyhedron object
        % 
        input_space

        % LTVSystem/input_dim
        % ====================================================================
        % 
        % Dimension of the input vector
        % 
        input_dim

        % LTVSystem/dist
        % ====================================================================
        % 
        % LTV system disturbance, RandomVector object
        % 
        dist

        % LTVSystem/dist_mat
        % ====================================================================
        % 
        % System disturbance matrix
        %   - MATLAB matrix for LTI system
        %   - Anonymous functions for LTV system, maps time index to matrix
        % 
        dist_mat

        % LTVSystem/dist_dim
        % ====================================================================
        % 
        % Dimension of disturbance vector
        % 
        dist_dim

    end

    properties (Access = private)
        % LTVSystem/sys_type
        % ====================================================================
        % 
        % Type of system; either LTI or LTV
        % 
        sys_type
    end

    methods
        function obj = LtvSystem(varargin)
        % Create a discrete-time LTI system object
        % ====================================================================
        %
        % Constructor method fot the LTV System class. Will create the 
        % LtvSystem object
        %
        % Usage:
        % ------
        % T = 0.5;
        % sys = LtvSystem('StateMatrix', [1, T; 0, 1], ...
        %                 'InputMatrix', [T^2/2;T], ...
        %                 'InputSpace', Polyhedron('lb', -1, 'ub', 1), ...
        %                 'DisturbanceMatrix', [T^2/2;T], ...
        %                 'Disturbance', Polyhedron('lb', -1, 'ub', 1));
        %
        % =====================================================================
        %
        % obj = LtvSystem(Name, Value)
        % 
        % Inputs:
        % -------
        %   ------------------------------------------------------------
        %   Name               | Value
        %   ------------------------------------------------------------
        %   StateMatrix        | Square numeric matrix
        %   InputMatrix        | (optional) Numeric matrix
        %   DisturbanceMatrix  | (optional) Numeric matrix
        %   InputSpace         | (optional) Polyhedron or empty
        %   Disturbance        | (optional) Polyhedron or RandomVector or empty
        % 
        % Outputs:
        % --------
        %   obj - LtvSystem object
        %
        % Notes:
        % ------
        % * 'InputMatrix' and 'InputSpace' need to be either defined together
        %   or neither of them.
        % * 'DisturbanceMatrix' and 'Disturbance' need to be either defined
        %   together or neither of them.
        % 
        % =====================================================================
        % 
        %   This function is part of the Stochastic Reachability Toolbox.
        %   License for the use of this function is given in
        %        https://sreachtools.github.io/license/
        % 

            inpar = inputParser();
            inpar.addParameter('StateMatrix', [], ...
                @(x) validateattributes(x, {'function_handle', 'numeric'}, ...
                    {'nonempty', 'square'}));
            inpar.addParameter('InputMatrix', [], ...
                @(x) validateattributes(x, {'function_handle', 'numeric'}, ...
                    {'nonempty'}));
            inpar.addParameter('InputSpace', Polyhedron(), ...
                @(x) isa(x, 'Polyhedron'));
            inpar.addParameter('Disturbance', [], ...
                @(x) validateattributes(x, {'RandomVector', 'Polyhedron'}, ...
                    {'nonempty'}));
            inpar.addParameter('DisturbanceMatrix', [], ...
                @(x) validateattributes(x, {'function_handle', 'numeric'}, ...
                    {'nonempty'}));

            try
                inpar.parse(varargin{:});
            catch err
                exc = SrtInvalidArgsError(['Invalid arguments provided to ', ...
                    'LtvSystem']);
                exc = exc.addCause(err);
                throwAsCaller(exc)
            end

            % Check that input matrix and space/disturbance matrix and space
            % are defined simultaneously
            if xor(any(strcmp(inpar.UsingDefaults, 'InputMatrix')), ...
                   any(strcmp(inpar.UsingDefaults, 'InputSpace')))

                exc = SrtInvalidArgsError(['InputMatrix and InputSpace ', ...
                    'must be simultaneously defined']);
                throw(exc)
            end

            if xor(any(strcmp(inpar.UsingDefaults, 'DisturbanceMatrix')), ...
                   any(strcmp(inpar.UsingDefaults, 'Disturbance')))

                exc = SrtInvalidArgsError(['Disturbance and ', ...
                    'DisturbanceMatrix must be simultaneously defined']);
                throw(exc)
            end

            % State matrix must be defined
            if any(strcmp(inpar.UsingDefaults, 'StateMatrix'))
                exc = SrtInvalidArgsError(['Invalid arguments provided to ', ...
                    'LtvSystem']);
                exc = exc.addCause(SrtInvalidArgsError(['StateMatrix must ', ...
                    'be provided to an LtvSystem']));
                throwAsCaller(exc)
            end

            % If all given matrices are numeric (empty is numeric) then system
            % is LTI
            if isnumeric(inpar.Results.StateMatrix) && ...
               isnumeric(inpar.Results.InputMatrix) && ...
               isnumeric(inpar.Results.DisturbanceMatrix)

                obj.sys_type = 'LTI';
            else
                obj.sys_type = 'LTV';
            end

            obj.input_space = inpar.Results.InputSpace;
            obj.input_dim = obj.input_space.Dim;

            obj.dist = inpar.Results.Disturbance;
            if any(strcmp(inpar.UsingDefaults, 'Disturbance'))
                obj.dist_dim = 0;
            else
                % Can be only Polyhedron or RandomVector due to input parsing
                if isa(obj.dist, 'RandomVector') 
                    obj.dist_dim = obj.dist.dim;
                else
                    obj.dist_dim = obj.dist.Dim;
                end
            end

            % Set the state matrix
            if isnumeric(inpar.Results.StateMatrix)
                if obj.islti()
                    obj.state_mat = inpar.Results.StateMatrix;
                else
                    obj.state_mat = @(t) inpar.Results.StateMatrix;
                end
            else
                obj.state_mat = inpar.Results.StateMatrix;
            end

            if obj.islti()
                A = obj.state_mat;
            else
                A = obj.state_mat(1);
            end

            % Set the state dimension size
            obj.state_dim = size(A, 1);

            % check to make sure that the obtained state matrices are square
            if size(A, 1) ~= size(A, 2)
                exc = SrtInvalidArgsError(['Invalid arguments provided to ', ...
                    'LtvSystem']);
                exc = exc.addCause(MException('', ...
                    'Obtained state matrix is not square'));
                throwAsCaller(exc)
            end


            %% Input matrix
            % If not specified then the input matrix should be a 
            % obj.state_dim x 0 matrix
            if any(strcmp(inpar.UsingDefaults, 'InputMatrix'))
                if obj.islti()
                    obj.input_mat = zeros(obj.state_dim, obj.input_dim);
                else
                    obj.input_mat = @(t) zeros(obj.state_dim, ...
                        obj.input_dim);
                end
            else
                if isnumeric(inpar.Results.InputMatrix)
                    if obj.islti()
                        obj.input_mat = inpar.Results.InputMatrix;
                    else
                        obj.input_mat = @(t) inpar.Results.InputMatrix;
                    end
                else
                    obj.input_mat = inpar.Results.InputMatrix;
                end
            end

            if obj.islti()
                B = obj.input_mat;
            else
                B = obj.input_mat(1);
            end

            % check consistency of the input matrix and state dimension, i.e.
            % that the input matrix has the same number of rows as the state
            % dimension
            if size(B, 1) ~= obj.state_dim
                exc = SrtInvalidArgsError(['Invalid arguments provided to ', ...
                    'LtvSystem']);
                exc = exc.addCause(MException('', ['Obtained input matrix ', ...
                    'is not consistent with the dimension of the ', ...
                    'state space, i.e. size(obj.input_mat(1), 1) ~= ', ...
                    'obj.state_dim']));
                throwAsCaller(exc)
            end

            % check consistency of the input matrix and input space dimension, 
            % i.e. that the input matrix has the same number of columns as the
            % dimension of the input space
            if size(B, 2) ~= obj.input_dim
                exc = SrtInvalidArgsError( ...
                    'Invalid arguments provided to LtvSystem');
                exc = exc.addCause(MException('', ['Obtained input matrix ', ...
                    'is not consistent with the dimension of the ', ...
                    'input space, i.e. size(obj.input_mat(1), 2) ~= ', ...
                    'obj.input.dim']));
                throwAsCaller(exc)
            end

            % Disturbance matrix
            % If not specified then the disturbance matrix should be a 
            % obj.state_dim x 0 matrix
            if any(strcmp(inpar.UsingDefaults, 'DisturbanceMatrix'))
                if obj.islti()
                    obj.dist_mat = zeros(obj.state_dim, obj.dist_dim);
                else
                    obj.dist_mat = @(t) zeros(obj.state_dim, ...
                        obj.dist_dim);
                end
            else
                if isnumeric(inpar.Results.DisturbanceMatrix)
                    if obj.islti()
                        obj.dist_mat = inpar.Results.DisturbanceMatrix;
                    else
                        obj.dist_mat = @(t) inpar.Results.DisturbanceMatrix;
                    end
                else
                    obj.dist_mat = inpar.Results.DisturbanceMatrix;
                end
            end
            
            if obj.islti()
                F = obj.dist_mat;
            else
                F = obj.dist_mat(1);
            end

            % check consistency of the input matrix and state dimension, i.e.
            % that the input matrix has the same number of rows as the state
            % dimension
            if size(F, 1) ~= obj.state_dim
                exc = SrtInvalidArgsError( ...
                    'Invalid arguments provided to LtvSystem');
                exc = exc.addCause(MException('', ['Obtained disturbance ', ...
                    'matrix is not consistent with the dimension of the ', ...
                    'state space, i.e. size(obj.dist_mat(1), 1) ~= ', ...
                    'obj.state_dim']));
                throwAsCaller(exc)
            end

            % check consistency of the disturbance matrix and input space dimension, 
            % i.e. that the disturbance matrix has the same number of columns as the
            % dimension of the input space
            if size(F, 2) ~= obj.dist_dim
                exc = SrtInvalidArgsError( ...
                    'Invalid arguments provided to LtvSystem');
                exc = exc.addCause(MException('', ['Obtained disturbance ', ...
                    'matrix is not consistent with the dimension of the ', ...
                    'disturbance space, i.e. size(obj.dist_mat(1), 2) ~= ', ...
                    'obj.dist_dim']));
                throwAsCaller(exc)
            end
        end

        function varargout = subsref(obj, s)
        % Overload of MATLAB internal subsref
        % ====================================================================
        %
        % Overloaded method of MATLAB's internal subsref. Overloading allows
        % for calling the functions (if existing) for the time-varying system
        % matrices
        %
        % Usage: Overload of internal method
        %
        % =====================================================================
        %
        % obj.property
        % obj.method(args)
        % obj.state_mat(t)
        % 
        % Inputs: None
        % 
        % Outputs:
        % --------
        %   v - Value from the subsref call
        %
        % =====================================================================
        % 
        %   This function is part of the Stochastic Reachability Toolbox.
        %   License for the use of this function is given in
        %        https://sreachtools.github.io/license/
        % 

            if  length(s) == 2 && ...
                strcmp(s(1).type, '.') && ...
                any(strcmp(s(1).subs, {'state_mat', 'input_mat', ...
                    'dist_mat'})) && ...
                strcmp(s(2).type, '()')

                % access matrix
                % check the paranthetical subscript indices
                if length(s(2).subs) > 1
                    exc = SrtInvalidArgsError( ...
                        ['Too many input arguments']);
                    throw(exc);
                end

                if ~isscalar(s(2).subs{1})
                    exc = SrtInvalidArgsError( ...
                        ['Time value must be scalar numeric input']);
                    throw(exc);
                end

                if obj.islti()
                    [varargout{1:nargout}] = obj.(s(1).subs);
                else
                    [varargout{1:nargout}] = obj.(s(1).subs)(s(2).subs{1});
                end
            else
                [varargout{1:nargout}] = builtin('subsref', obj, s);
            end
        end

        function disp(obj, varargin)
        % Overload of MATLAB internal disp
        % ====================================================================
        %
        % Overloaded method of MATLAB's internal disp. 
        %
        % Usage: Overload of internal method
        %
        % =====================================================================
        %
        % disp(obj, Name, Value)
        % 
        % Inputs:
        % -------
        %   obj - LtvSystem object
        %   ------------------------------------------------------------
        %   Name           | Value
        %   ------------------------------------------------------------
        %   verbose        | true or false
        % 
        % Outputs: None
        % 
        % Notes:
        % ------
        % * disp function for this class was inspired from MPT3
        %   (http://people.ee.ethz.ch/~mpt/3/)
        %
        % =====================================================================
        % 
        %   This function is part of the Stochastic Reachability Toolbox.
        %   License for the use of this function is given in
        %        https://sreachtools.github.io/license/
        % 

            inpar = inputParser();
            inpar.addParameter('verbose', false, ...
                @(x) assert(islogical(x), 'verbose input must be a logical'));

            try
                inpar.parse(varargin{:})
            catch err
                exc = MException('SReachTools:innvalidArgs', ...
                    'Invalid arguments to LtvSystem/disp');
                exc = exc.addCause(err);
                throwAsCaller(exc)
            end

            fprintf(['Linear time varying system with %d states, ', ...
                '%d inputs, and %d disturbances.\n'], obj.state_dim, ...
                obj.input_dim, obj.dist_dim);
            if inpar.Results.verbose
                fprintf('\nState matrix function -> %s\n\n', ...
                    func2str(obj.state_mat))
                if obj.input_dim > 0
                    disp('Input space')
                    disp(obj.input_space)
                    fprintf('Input matrix function -> %s\n\n', ...
                        func2str(obj.input_mat))
                else
                    fprintf('System has no input\n\n')
                end
                if obj.dist_dim > 0
                    if isa(obj.dist, 'Polyhedron')
                        fprintf('Disturbance is of type: %s\n', class(obj.dist))
                    else
                        fprintf('Disturbance is a ')
                        disp(obj.dist)
                    end
                    fprintf('Disturbance matrix function -> %s\n\n', ...
                        func2str(obj.dist_mat))
                else
                    fprintf('System has no disturbance\n\n')
                end
            end
        end

        function yn = islti(obj)
        % Get boolean result if system is LTI
        % ====================================================================
        %
        % Get boolean result if system is LTI. Considered LTI if state_mat, 
        % input_mat, and dist_mat are all matrices
        %
        % Usage:
        % ------
        % T = 0.5;
        % sys = LtvSystem('StateMatrix', [1, T; 0, 1], ...
        %                 'InputMatrix', [T^2/2;T], ...
        %                 'InputSpace', Polyhedron('lb', -1, 'ub', 1), ...
        %                 'DisturbanceMatrix', [T^2/2;T], ...
        %                 'Disturbance', Polyhedron('lb', -1, 'ub', 1));
        % if sys.islti()
        %   disp('System is LTI')
        % else
        %   disp('System is not LTI')
        % end
        % 
        % =====================================================================
        %
        % yn = obj.islti()
        % 
        % Inputs: None
        % Outputs:
        % --------
        %   yn - Logical value of 1 if system is LTI
        % 
        % =====================================================================
        % 
        %   This function is part of the Stochastic Reachability Toolbox.
        %   License for the use of this function is given in
        %        https://sreachtools.github.io/license/
        % 

            yn = strcmp(obj.sys_type, 'LTI');
        end

        function yn = isltv(obj)
        %  Get boolean result if system is LTI
        % ====================================================================
        %
        % Get boolean result if system is LTI. Considered LTI if state_mat, 
        % input_mat, and dist_mat are all matrices
        %
        % Usage:
        % ------
        % T = 0.5;
        % sys = LtvSystem('StateMatrix', [1, T; 0, 1], ...
        %                 'InputMatrix', [T^2/2;T], ...
        %                 'InputSpace', Polyhedron('lb', -1, 'ub', 1), ...
        %                 'DisturbanceMatrix', [T^2/2;T], ...
        %                 'Disturbance', Polyhedron('lb', -1, 'ub', 1));
        % if sys.islti()
        %   disp('System is LTI')
        % else
        %   disp('System is not LTI')
        % end
        % 
        % =====================================================================
        %
        % yn = obj.islti()
        % 
        % Inputs: None
        % Outputs:
        % --------
        %   yn - Logical value of 1 if system is LTI
        % 
        % =====================================================================
        % 
        %   This function is part of the Stochastic Reachability Toolbox.
        %   License for the use of this function is given in
        %        https://sreachtools.github.io/license/
        % 

            yn = strcmp(obj.sys_type, 'LTV');
        end
        
        % Methods that have been defined externally
        [Z, H, G] = getConcatMats(sys, time_horizon);
        [concat_input_space_A, concat_input_space_b] = ...
            getConcatInputSpace(sys, time_horizon);
    end
end
