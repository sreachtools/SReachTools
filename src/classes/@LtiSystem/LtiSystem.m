classdef LtiSystem
% SReachTools/LtiSystem: Create a discrete-time LTI system object
% ============================================================================
%
% Defines a discrete-time LTI system that is:
%     - control-free and disturbance-free, or
%     - controlled but disturbance-free, or
%     - perturbed but control-free, or
%     - controlled and perturbed system.
%
% Perturbation can be either:
%     - a bounded uncertainity with no stochastic information
%     - a StochasticDisturbance object
%
%  Usage:
%  ------
%  % Define a double integrator system:
%
%  T = 0.5;
%  sys = LtiSystem('StateMatrix', [1, T; 0, 1], ...
%                  'InputMatrix', [T^2/2;T], ...
%                  'InputSpace', Polyhedron('lb', -1, 'ub', 1), ...
%                  'DisturbanceMatrix', [T^2/2;T], ...
%                  'Disturbance', Polyhedron('lb', -1, 'ub', 1));
%   
% ============================================================================
%
% LTISYSTEM Properties:
% ---------------------
%   state_matrix          - State matrix (Square matrix)
%   input_matrix          - Input matrix
%   input_space           - Input space (Polyhedron)
%   disturbance_matrix    - Disturbance matrix  
%   disturbance           - Disturbance object (empty/Polyhedron/StochasticDisturbance)     
%   state_dimension       - State dimension (scalar)   
%   input_dimension       - Input dimension (scalar)  
%   disturbance_dimension - Disturbance dimension (scalar)
% 
% LTISYSTEM Methods:
% ------------------
%   LtiSystem/LtiSystem   - Constructor
%   getConcatInputSpace   - Get concatenated input space
%   getConcatMats         - Get concatenated input and state matrices
% 
%
%
% Notes:
% -----
% * EXTERNAL DEPENDENCY: Uses MPT3 to define input,robust disturbance space
% =============================================================================
%
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/abyvinod/SReachTools/blob/master/LICENSE
%
%

    properties (SetAccess = immutable)
        state_matrix          = []
        state_dimension       = []
        input_matrix          = []
        input_space           = []
        input_dimension       = []
        disturbance_matrix    = []
        disturbance           = []
        disturbance_dimension = []
    end
    
    methods
        function obj = LtiSystem(varargin)
        % SReachTools/LtiSystem/LtiSystem: Create a discrete-time LTI system 
        % object
        % ====================================================================
        %
        % Constructor method fot the LTI System class. Will create the 
        % LtiSystem object
        %
        % Usage:
        % ------
        % T = 0.5;
        % sys = LTISYSTEM('StateMatrix', [1, T; 0, 1], ...
        %                 'InputMatrix', [T^2/2;T], ...
        %                 'InputSpace', Polyhedron('lb', -1, 'ub', 1), ...
        %                 'DisturbanceMatrix', [T^2/2;T], ...
        %                 'Disturbance', Polyhedron('lb', -1, 'ub', 1));
        %
        % =====================================================================
        %
        % obj = LTISYSTEM(Name, Value)
        % 
        % Inputs:
        % -------
        %   ------------------------------------------------------------
        %   Name               | Value
        %   ------------------------------------------------------------
        %   StateMatrix        | Square numeric matrix
        %   InputMatrix        | (optional) Numeric matrix
        %   DisturbanceMatrix  | (optional) Numeric matrix
        %   InputSpace         | (optional) Polyhedron
        %   Disturbance        | (optional) Polyhedron or 
        %                      |            StochasticDisturbance
        % 
        % Outputs:
        % --------
        %   obj - LtiSystem object
        %
        % Notes:
        %   - 'InputMatrix' and 'InputSpace' need to be either defined together or
        %     neither of them.
        %   - 'DisturbanceMatrix' and 'Disturbance' need to be either defined together 
        %     or neither of them.
        % 
        % =====================================================================
        % 
        %   This function is part of the Stochastic Reachability Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/abyvinod/SReachTools/blob/master/LICENSE
        % 

            % Input arguments must be name-value pairs (hence even in count)
            assert(mod(length(varargin), 2) == 0, ...
                   'SReachTools:invalidArgs', ...  
                   'Arguments must be given as name-value pairs');
            
            var_index = 1;
            while var_index < length(varargin)
                % Loop through the input arguments (Name followed by value)
                switch(lower(varargin{var_index}))
                    case 'statematrix'
                        validateattributes(varargin{var_index+1}, ...
                                          {'numeric'}, ...
                                          {'nonempty'});
                        obj.state_matrix = varargin{var_index+1};
                    case 'inputmatrix'
                        validateattributes(varargin{var_index+1}, ...
                                           {'numeric'}, ...
                                           {'nonempty'});
                        obj.input_matrix = varargin{var_index+1};
                    case 'disturbancematrix'
                        validateattributes(varargin{var_index+1}, ...
                                           {'numeric'}, ...
                                           {'nonempty'});
                        obj.disturbance_matrix = varargin{var_index+1};
                    case 'inputspace'
                        % Currently InputSpace has to be a Polyhedron (MPT)
                        assert(exist('mpt_init','file')==2, ...
                               'SReachTools:setup_error', ...
                               ['This function uses MPT3. Please get it ', ...
                                'from http://control.ee.ethz.ch/~mpt/3/.']);
                        validateattributes(varargin{var_index+1}, ...
                                           {'Polyhedron'}, ...
                                           {'nonempty'});
                        obj.input_space = varargin{var_index+1};
                    case 'disturbance'
                        % Currently Disturbance has to be a Polyhedron (MPT) or
                        % an object of StochasticDisturbance
                        validateattributes(varargin{var_index+1}, ...
                                           {'Polyhedron', ...
                                            'StochasticDisturbance'}, ...
                                           {'nonempty'});
                        obj.disturbance = varargin{var_index+1};
                    otherwise
                        % Raise exception for any unhandled argument
                        error('SReachTools:invalidArgs', ...
                              'Unhandled argument given');
                end
                var_index = var_index + 2;
            end

            % Defining state dimension after ensuring a non-empty square state
            % matrix
            assert(~isempty(obj.state_matrix), ...
                   'SReachTools:invalidArgs', ...  
                   'State matrix can not be empty');
            assert(size(obj.state_matrix,2) == size(obj.state_matrix,1), ...
                   'SReachTools:invalidArgs', ...
                   'State matrix is not square');
            obj.state_dimension = size(obj.state_matrix,2);

            % Setting default values for properties that were not specified
            % No values provided to matrices => zero
            if isempty(obj.input_matrix)
                obj.input_matrix = zeros(obj.state_dimension,1);
            end
            if isempty(obj.disturbance_matrix)
                obj.disturbance_matrix = zeros(obj.state_dimension,1);
            end
            % No values provided to 'spaces' => empty polyhedra
            if isempty(obj.input_space)
                obj.input_space = Polyhedron();
            end 
            if isempty(obj.disturbance)
                obj.disturbance = Polyhedron();
            end 

            % Updating dimension information
            obj.input_dimension = obj.input_space.Dim;
            if strcmp(class(obj.disturbance), 'Polyhedron')
                obj.disturbance_dimension = obj.disturbance.Dim;
            elseif strcmp(class(obj.disturbance), 'StochasticDisturbance')
                obj.disturbance_dimension = obj.disturbance.dimension;
            else
                error('SReachTools:internal', ...
                      'Unsupported disturbance provided');
            end

            % Sanity checks
            checkSystemProperties(obj);
        end

        function disp(obj)
        % SReachTools/LtiSystem/disp: Display information about LtiSystem object
        % =====================================================================
        %
        % Overloading of MATLAB's built-in display for object to create
        % prettier and more concise output
        %
        % Usage
        % -----
        % % needs variable 'ltisys' which is an LtiSystem object
        % disp(ltisys);
        % ltisys
        % ltisys.disp();
        % 
        % =====================================================================
        % 
        % disp(obj)
        % 
        % Inputs: None
        % Outputs: None
        %
        % Notes:
        %   - disp function for this class was inspired from MPT3
        %     (http://people.ee.ethz.ch/~mpt/3/)
        % 
        % =====================================================================
        % 
        %   This function is part of the Stochastic Reachability Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/abyvinod/SReachTools/blob/master/LICENSE
        %           

            plural = @(s, n) [num2str(n) ' ' s repmat('s', 1, double(n~=1))];
            disp(sprintf('LTI System with %s, %s, %s', ...
				plural('state', obj.state_dimension), ...
				plural('input', obj.input_dimension), ...
				plural('disturbance', obj.disturbance_dimension)));
        end
    end

    
    methods (Hidden)
        function checkSystemProperties(obj)
        % SReachTools/LtiSystem/checkSystemProperties: Check the properties of
        % and LtiSystem object
        % =====================================================================
        %
        % Class (hidden) method to verify the correct definition of system 
        % properties for an LtiSystem object
        %
        % Usage:
        % ------
        % % needs variable 'ltisys' which is an LtiSystem object
        % ltisys.checkSystemProperties();
        % 
        % ====================================================================
        % 
        % checkSystemProperties(obj)
        % 
        % Inputs: None
        % Outputs: None
        %
        % =====================================================================
        % 
        %   This function is part of the Stochastic Reachability Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/abyvinod/SReachTools/blob/master/LICENSE
        %   
            
            % Developer note: Joseph Gleason - 30 March 2018
            % This method should probably be changed to a private method

            % Sanity check --- Have all properties been initialized?
            props = properties(obj);
            for i = 1:length(props)
                if isempty(obj.(props{i}))
                    error('SReachTools:internal', ...
                          sprintf('Property %s not init by the constructor', ...
                                  props{i}));
                end
            end
            % Sanity check --- Input matrix of correct size
            assert(obj.state_dimension == size(obj.input_matrix,1), ...
                   'SReachTools:invalidArgs', ...
                   'Input matrix does not have correct row numbers');
            if obj.input_dimension > 0
                assert(obj.input_dimension == size(obj.input_matrix,2), ...
                       'SReachTools:invalidArgs', ...
                       'Input matrix does not have correct column numbers');
            else
                assert(obj.input_dimension + 1 == size(obj.input_matrix,2), ...
                       'SReachTools:invalidArgs', ...
                       'Empty input space: But non-column input matrix');
            end
            % Sanity check --- Disturbance matrix of correct size
            assert(obj.state_dimension == size(obj.disturbance_matrix,1), ...
                   'SReachTools:invalidArgs', ...
                   'Disturbance matrix does not have correct row numbers');
            if obj.disturbance_dimension > 0
                % Empty polyhedron has dimension 0, but matrix set to [0]
                assert(obj.disturbance_dimension == ...
                       size(obj.disturbance_matrix,2), ...
                       'SReachTools:invalidArgs', ...
                       ['Disturbance matrix does not have correct column', ... 
                       ' numbers']);
            else
                % Empty polyhedron has dimension 0, but matrix set to [0]
                assert(obj.disturbance_dimension + 1 == ...
                       size(obj.disturbance_matrix,2), ...
                       'SReachTools:invalidArgs', ...
                       'Empty disturbance: But non-column disturbance matrix');
            end
        end
    end
end
