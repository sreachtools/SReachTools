classdef LtvSystem 

    properties (SetAccess = protected)
        state_mat
        state_dim

        input_mat
        input_space
        input_dim

        disturbance
        disturbance_mat
        dist_dim
        type
    end

    methods
        function obj = LtvSystem(varargin)

            inpar = inputParser();
            inpar.addParameter('StateMatrix', [], ...
                @(x) validateattributes(x, {'function_handle', 'numeric'}, ...
                    {'nonempty', 'square'}));
            inpar.addParameter('InputMatrix', [], ...
                @(x) validateattributes(x, {'function_handle', 'numeric'}, ...
                    {'nonempty', 'square'}));
            inpar.addParameter('InputSpace', Polyhedron(), ...
                @(x) isa(x, 'Polyhedron'));
            inpar.addParameter('Disturbance', [], ...
                @(x) validateattributes(x, {'RandomVector', 'Polyhedron'}, ...
                    {'nonempty'}));
            inpar.addParameter('DisturbanceMatrix', [], ...
                @(x) validateattributes(x, {'function_handle', 'numeric'}, ...
                    {'nonempty', 'square'}));

            try
                inpar.parse(varargin{:});
            catch err
                exc = MException('SReachTools:invalidArgs', ...
                    'Invalid arguments provided to LtvSystem');
                exc = exc.addCause(err);
                throwAsCaller(exc)
            end

            % Check that input matrix and space/disturbance matrix and space
            % are defined simultaneously
            if xor(any(strcmp(inpar.UsingDefaults, 'InputMatrix')), ...
                   any(strcmp(inpar.UsingDefaults, 'InputSpace')))

                exc = MException('SReachTools:invalidArgs', ...
                    ['InputMatrix and InputSpace must be simultaneously ', ...
                     'defined']);
                throw(exc)
            end

            if xor(any(strcmp(inpar.UsingDefaults, 'DisturbanceMatrix')), ...
                   any(strcmp(inpar.UsingDefaults, 'Disturbance')))

                exc = MException('SReachTools:invalidArgs', ...
                    ['Disturbance and DisturbanceMatrix must be simultaneously ', ...
                     'defined']);
                throw(exc)
            end

            % State matrix must be defined
            if any(strcmp(inpar.UsingDefaults, 'StateMatrix'))
                exc = MException('SReachTools:invalidArgs', ...
                    'Invalid arguments provided to LtvSystem');
                exc = exc.addCause(MException('SReachTools:invalidArgs', ...
                    ['StateMatrix must be provided to an LtvSystem']));
                throwAsCaller(exc)
            end

            % If all given matrices are numeric (empty is numeric) then system
            % is LTI
            if isnumeric(inpar.Results.StateMatrix) && ...
               isnumeric(inpar.Results.InputMatrix) && ...
               isnumeric(inpar.Results.DisturbanceMatrix)

                obj.type = 'LTI';
            else
                obj.type = 'LTV';
            end

            obj.input_space = inpar.Results.InputSpace;
            obj.input_dim = obj.input_space.Dim;

            obj.disturbance = inpar.Results.Disturbance;
            if any(strcmp(inpar.UsingDefaults, 'Disturbance'))
                obj.dist_dim = 0;
            else
                obj.dist_dim = obj.disturbance.dimension;
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
                exc = MException('SReachTools:invalidArgs', ...
                    'Invalid arguments provided to LtvSystem');
                exc = exc.addCause(MException('', 'Obtained state matrix is not square'));
                throwAsCaller(exc)
            end


            % Input matrix
            % If not specified then the input matrix should be a 
            % obj.state_dim x 0 matrix
            if isnumeric(inpar.Results.InputMatrix)
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
                exc = MException('SReachTools:invalidArgs', ...
                    'Invalid arguments provided to LtvSystem');
                exc = exc.addCause(MException('', ['Obtained input matrix is not ', ...
                    'consistent with the dimension of the state space, ', ...
                    'i.e. size(obj.input_mat(1), 1) ~= obj.state_dim']));
                throwAsCaller(exc)
            end

            % check consistency of the input matrix and input space dimension, 
            % i.e. that the input matrix has the same number of columns as the
            % dimension of the input space
            if size(B, 2) ~= obj.input_dim
                exc = MException('SReachTools:invalidArgs', ...
                    'Invalid arguments provided to LtvSystem');
                exc = exc.addCause(MException('', ['Obtained input matrix is not ', ...
                    'consistent with the dimension of the input space, ', ...
                    'i.e. size(obj.input_mat(1), 2) ~= obj.input.dim']));
                throwAsCaller(exc)
            end

            % Input matrix
            % If not specified then the input matrix should be a 
            % obj.state_dim x 0 matrix
            if isnumeric(inpar.Results.DisturbanceMatrix)
                if any(strcmp(inpar.UsingDefaults, 'DisturbanceMatrix'))
                    if obj.islti()
                        obj.disturbance_mat = zeros(obj.state_dim, obj.input_dim);
                    else
                        obj.disturbance_mat = @(t) zeros(obj.state_dim, ...
                            obj.input_dim);
                    end
                else
                    if isnumeric(inpar.Results.DisturbanceMatrix)
                        if obj.islti()
                            obj.disturbance_mat = inpar.Results.InputMatrix;
                        else
                            obj.disturbance_mat = @(t) inpar.Results.InputMatrix;
                        end
                    else
                        obj.disturbance_mat = inpar.Results.DisturbanceMatrix;
                    end
                end
            end

            if obj.islti()
                F = obj.disturbance_mat;
            else
                F = obj.disturbance_mat(1);
            end

            % check consistency of the input matrix and state dimension, i.e.
            % that the input matrix has the same number of rows as the state
            % dimension
            if size(F, 1) ~= obj.state_dim
                exc = MException('SReachTools:invalidArgs', ...
                    'Invalid arguments provided to LtvSystem');
                exc = exc.addCause(MException('', ['Obtained disturbance matrix is not ', ...
                    'consistent with the dimension of the state space, ', ...
                    'i.e. size(obj.disturbance_mat(1), 1) ~= obj.state_dim']));
                throwAsCaller(exc)
            end

            % check consistency of the disturbance matrix and input space dimension, 
            % i.e. that the disturbance matrix has the same number of columns as the
            % dimension of the input space
            if size(F, 2) ~= obj.dist_dim
                exc = MException('SReachTools:invalidArgs', ...
                    'Invalid arguments provided to LtvSystem');
                exc = exc.addCause(MException('', ['Obtained disturbance matrix is not ', ...
                    'consistent with the dimension of the disturbance space, ', ...
                    'i.e. size(obj.disturbance_mat(1), 2) ~= obj.input.dim']));
                throwAsCaller(exc)
            end
        end

        function v = subsref(obj, s)
            if  length(s) == 2 && ...
                strcmp(s(1).type, '.') && ...
                any(strcmp(s(1).subs, {'state_mat', 'input_mat', ...
                    'disturbance_mat'})) && ...
                strcmp(s(2).type, '()')

                % access matrix
                % check the paranthetical subscript indices
                if length(s(2).subs) > 1
                    exc = MException('SReachTools:invalidArgs', ...
                        ['Too many input arguments']);
                    throw(exc);
                end

                if ~isscalar(s(2).subs{1})
                    exc = MException('SReachTools:invalidArgs', ...
                        ['Time value must be scalar numeric input']);
                    throw(exc);
                end

                if obj.islti()
                    v = obj.(s(1).subs);
                else
                    v = obj.(s(1).subs)(s(2).subs{1});
                end
            else
                v = builtin('subsref', obj, s);
            end
        end

        function disp(obj, varargin)

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

            if inpar.Results.verbose
                fprintf(['Linear Time Varying System with:\n\n', ...
                    '    %d states\n', ...
                    '    %d inputs\n', ...
                    '    state matrix function -> %s\n', ...
                    '    input matrix function -> %s\n\n'], ...
                    obj.state_dim, ...
                    obj.input_dim, ...
                    length(obj.params{1}), ...
                    func2str(obj.state_mat), ...
                    func2str(obj.input_mat));
            else
                fprintf(['Linear Time Varying System with %d states, ', ...
                    'and %d inputs\n'], obj.state_dim, ...
                    obj.input_dim);
            end
        end

        function ans = islti(obj)
            ans = strcmp(obj.type, 'LTI');
        end

        function ans = isltv(obj)
            ans = strcmp(obj.type, 'LTV');
        end
    end

    methods (Static, Access = private)
        function validateParameters(params)
            if ~isa(params, 'cell') || ~isvector(params)
                error('SReachTools:invalidArgs', ['Parameters must be ', ...
                    'a cell array of cell arrays']);
            end

            for lv = 1:length(params)
                if ~isa(params{lv}, 'cell')
                    error('SReachTools:invalidArgs', ['Parameters must be ', ...
                        'a cell array of cell arrays']);
                end
            end
        end
    end

end