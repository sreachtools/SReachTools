classdef RandomVector
% Create a random vector object
% ==========================================================================
%
% Defines a random vector with a standard probability density function (pdf)
%
% We support the following pdfs:
%     1. Gaussian
%
% Usage:
% ------
%
% % Define a Gaussian random variable of mean 0 and standard deviation 2:
% GaussianRV = RandomVector('Gaussian', ...
%                           0, ...
%                           2^2);
%
% % Define a Gaussian random vector of mean [0;2] and covariance matrix 
% % eye(2):
% GaussianRV = RandomVector('Gaussian', ...
%                           [0;2], ...
%                           eye(2));
%   
% ==========================================================================
%
% RandomVector Properties:
% ------------------------
%   type       - Random vector type (string)
%   parameters - System parameters (struct)
%   dim        - Random vector dimension (scalar)
%   pdf        - Probability density function (function handle)
%
% RandomVector Methods:
% ---------------------
%   RandomVector/RandomVector - Class constructor
% 
% Notes:
% ------
% * MATLAB DEPENDENCY: Uses MATLAB's Statistics and Machine Learning Toolbox.
%                      Needs mvnpdf
% * Currently only supports Gaussian random vectors
% * Requires the mean and the covariance matrices to be non-empty column
%   vector and a symmetric matrix respectively
% * The anonymous function used for the definition of obj.pdf transposes the
%   accepted column vector for using mvnpdf.
% * RandomVector.pdf takes in arguments of the form N_points x
%   random_vector_dim
% 
% =========================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 
% 

    properties
        % RandomVector/type
        % ==================================================================
        % 
        % String indicator of type of random vector
        %
        % Acceptable types:
        %   Gaussian
        % 
        type

        % RandomVector/parameters
        % ==================================================================
        % 
        % Struct containing random vector parameter information; will vary for
        % different random vector types. Currently only Gaussian random vector
        % are supported.
        %
        % Gaussian type:
        %   parameters.mean       - Mean vector (p x 1)
        %   parameters.covariance - Covariance matrix (p x p)
        % 
        parameters

        % RandomVector/dim
        % ==================================================================
        % 
        % Dimension of the random vector
        % 
        dim

        % RandomVector/pdf
        % ==================================================================
        % 
        % Probability density function of random vector
        % 
        pdf
    end
    methods
        function obj = RandomVector(rv_type, varargin)
        %  Constructor for random vector class
        % ====================================================================
        %
        % Inputs:
        % -------
        %   rv_type  - Character array of random vector type; currently
        %              acceptable types:
        %                  'Gaussian'
        %   varargin - Arguments with the properties of the random vector. For
        %              Gaussian random vector:
        %                   mu
        %                   sigma
        %
        % Outputs:
        % --------
        %   obj - Random vector object
        %
        % =====================================================================
        % 
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        % 
            
            % Check if the random vector type is a string 
            validateattributes(rv_type, {'char'}, {'nonempty'});
            % Update the random vector type
            obj.type = rv_type; 

            switch(lower(obj.type))
                case 'gaussian'
                    if length(varargin) ~= 2 
                        throwAsCaller(SrtInvalidArgsError(['Gaussian random vector needs the mean vector ', ...
                           'and covariance matrix']));
                    end

                    % Check if mean is a nonempty numeric column vector
                    validateattributes(varargin{1}, {'numeric'}, {'column'});
                    % Check if covariance matrix is a nonempty square matrix
                    validateattributes(varargin{2}, {'numeric'}, ...
                                       {'nonempty','square'});
                                   
                    % set parameters
                    obj.parameters.mean = varargin{1};
                    obj.parameters.covariance = varargin{2};
                    
                    % TODO: Throw warning as in ellipsoid
                    if ~issymmetric(obj.parameters.covariance)
                        % Compute the symmetric component of it
                        symm_cov_matrix = (obj.parameters.covariance +...
                            obj.parameters.covariance')/2;
                        % Max error element-wise
                        max_err = max(max(abs(obj.parameters.covariance -...
                            symm_cov_matrix)));
                        if max_err > eps
                            warning('SReachTools:runtime',sprintf(...
                                ['Non-symmetric covariance matrix made ',...
                                 'symmetric (max element-wise error: ',...
                                 '%1.3e)!'], max_err));
                        end
                        obj.parameters.covariance = symm_cov_matrix;
                    end
                    
                    % Check if the mean and covariance are of correct dimensions
                    if size(obj.parameters.mean, 1) ~=...
                            size(obj.parameters.covariance, 1)
                        throwAsCaller(SrtInvalidArgsError(['Mean and ',...
                            'covariance matrix have different dimensions']));
                    end
                        
                    % For some reason, -eps alone is not enough?
                    min_eig_val = min(eig(obj.parameters.covariance));
                    if  min_eig_val < -2*eps
                        throwAsCaller(SrtInvalidArgsError(['Covariance ',...
                            'matrix can not have negative eigenvalues']));
                    elseif min_eig_val <= eps
                        warning('SReachTools:runtime',['Creating a',...
                            ' Gaussian which might have a deterministic ',...
                            'component']);
                    end
                                        
                    % Update the dimension
                    obj.dim = size(obj.parameters.mean, 1);
                    
                    % Define an anonymous function using MATLAB's built-in pdf 
                    % Transpose the mean for mvnpdf
                    obj.pdf = @(x) mvnpdf(x, ...
                                          obj.parameters.mean', ...
                                          obj.parameters.covariance);
                otherwise
                    exc = SrtInternalError(['Unsupported random ', ...
                        'vector type']);
                    throw(exc);
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
        %      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        %
            
            fprintf('%d-dimensional %s random vector\n', obj.dim, obj.type);
        end
        
        function newobj=mtimes(obj, F)
        % Override of MATLAB multiplication command
        % ====================================================================
        % 
        % Inputs:
        % -------
        %   obj - RandomVector object
        %   F   - Linear transformation matrix for multiplication
        %
        % Outputs:
        % --------
        %   newobj - RandomVector object (F*obj)
        %
        % Notes:
        % ------
        % * Requires Gaussian random vector assumption
        % ====================================================================
        % 
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        %
            
            switch [class(obj), class(F)]
                case ['RandomVector','double']
                    % All ok
                case ['double', 'RandomVector']
                    % Need to switch the arguments
                    Ftemp = obj;
                    obj = F;
                    F = Ftemp;
                otherwise
                    throwAsCaller(SrtInvalidArgsError(sprintf(['Operation *',...
                       ' not defined between *%s, %s'], class(obj), class(F))));
            end
            switch obj.type
                case 'Gaussian'
                    if size(F, 2) ~= obj.dim
                        throwAsCaller(SrtInvalidArgsError(['Mismatch ',...
                            'between dimensions']));
                    end
                    newobj=RandomVector('Gaussian', F*obj.parameters.mean,...
                        F*obj.parameters.covariance*F');
                otherwise
                    throwAsCaller(SrtInvalidArgsError(['Multiplication is ',...
                        'supported only for Gaussian random vectors']));
            end
        end
        
        function newobj=concat(obj, time_horizon)
        % Create a concatenated random vector of length time_horizon
        % ====================================================================
        % 
        % Inputs:
        % -------
        %   obj           - RandomVector object (typically disturbance w_k)
        %   time_horizon  - The number of time steps of interest N
        %
        % Outputs:
        % --------
        %   newobj        - RandomVector object (for disturbance, it can be used 
        %                   as W = [w_0^\top w_1^\top ... w_{N-1}^\top])
        %
        % Notes:
        % ------
        % * Requires Gaussian random vector assumption
        %
        % ====================================================================
        % 
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        %
            
            % Check if mean is a nonempty numeric column vector
            validateattributes(time_horizon, {'numeric'},...
                {'scalar','integer','>','0'});
            switch obj.type
                case 'Gaussian'
                    muW = repmat(obj.parameters.mean,time_horizon,1);
                    covW = kron(eye(time_horizon), obj.parameters.covariance);
                    newobj=RandomVector('Gaussian', muW, covW);
                otherwise
                    throwAsCaller(SrtInvalidArgsError(['Concatenation is ',...
                        'supported only for Gaussian random vectors']));
            end
        end
    end
end
