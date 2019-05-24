classdef LtiSystem < LtvSystem
% Create a discrete-time LTI system object
% ============================================================================
%
% Defines a discrete-time LTI system that is:
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
%   state_matrix - State matrix (Square matrix, state_dim x state_dim)
%   input_matrix - Input matrix (Matrix, state_dim x input_dim)
%   input_space  - Input space (empty / Polyhedron)
%   dist_matrix  - Disturbance matrix (Matrix, state_dim x dist_dim)
%   dist         - Disturbance object 
%                  (empty / Polyhedron / RandomVector)     
%   state_dim    - State dimension (scalar)   
%   input_dim    - Input dimension (scalar)  
%   dist_dim     - Disturbance dimension (scalar)
% 
% LTISYSTEM Methods:
% ------------------
%   LtiSystem/LtiSystem   - Constructor
%   getConcatInputSpace   - Get concatenated input space
%   getConcatMats         - Get concatenated state, input, and disturbance
%                           matrices
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
    
    methods
        function obj = LtiSystem(varargin)
        % Create a discrete-time LTI system object
        % ====================================================================
        %
        % Constructor method fot the LTI System class. Will create the 
        % LtiSystem object
        %
        % Usage:
        % ------
        % T = 0.5;
        % sys = LtiSystem('StateMatrix', [1, T; 0, 1], ...
        %                 'InputMatrix', [T^2/2;T], ...
        %                 'InputSpace', Polyhedron('lb', -1, 'ub', 1), ...
        %                 'DisturbanceMatrix', [T^2/2;T], ...
        %                 'Disturbance', Polyhedron('lb', -1, 'ub', 1));
        %
        % =====================================================================
        %
        % obj = LtiSystem(Name, Value)
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
        %   obj - LtiSystem object
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
                @(x) validateattributes(x, {'numeric'}, ...
                    {'nonempty', 'square'}));
            inpar.addParameter('InputMatrix', [], ...
                @(x) validateattributes(x, {'numeric'}, ...
                    {'nonempty'}));
            inpar.addParameter('DisturbanceMatrix', [], ...
                @(x) validateattributes(x, {'numeric'}, ...
                    {'nonempty'}));

            % Because InputSpace and Disturbance will be handled in the
            % LtvSystem superclass call just ignore what their inputs
            inpar.addParameter('InputSpace', Polyhedron(), @(x) true);
            inpar.addParameter('Disturbance', [], @(x) true);

            try
                inpar.parse(varargin{:});
            catch err
                exc = SrtInvalidArgsError(['Invalid arguments provided to ', ...
                    'LtiSystem']);
                exc = exc.addCause(err);
                throwAsCaller(exc)
            end

            obj@LtvSystem(varargin{:})            
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
        %   obj - LtiSystem object
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
                    'Invalid arguments to LtiSystem/disp');
                exc = exc.addCause(err);
                throwAsCaller(exc)
            end

            fprintf(['Linear time invariant system with %d states, ', ...
                '%d inputs, and %d disturbances.\n\n'], obj.state_dim, ...
                obj.input_dim, obj.dist_dim);
            if inpar.Results.verbose
                disp('State matrix')
                disp(obj.state_mat)
                if obj.input_dim > 0
                    disp('Input space')
                    disp(obj.input_space)
                    disp('Input matrix')
                    disp(obj.input_mat)
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
                    disp('Disturbance matrix')
                    disp(obj.dist_mat)
                else
                    fprintf('System has no disturbance\n\n')
                end
            end
        end

        % Methods inherited from LTV
        function [Z,G,H] = getConcatMats(sys, time_horizon)
            [Z,G,H] = getConcatMats@LtvSystem(sys, time_horizon);
        end

        function [concat_input_space_A, concat_input_space_b] =...
            getConcatInputSpace(sys, time_horizon)

            [concat_input_space_A, concat_input_space_b] =...
                getConcatInputSpace@LtvSystem(sys, time_horizon);
        end
    end
end
