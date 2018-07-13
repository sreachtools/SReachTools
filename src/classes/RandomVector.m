classdef RandomVector
% SReachTools/RandomVector: Create a random vector object
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
%                           0,
%                           2^2);
%
% % Define a Gaussian random vector of mean [0;2] and covariance matrix 
% % eye(2):
% GaussianRV = RandomVector('Gaussian', ...
%                           [0;2],
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
        dim
        pdf
    end
    methods
        function obj = RandomVector(rv_type, varargin)
        % SReachTools/RandomVector: Constructor for random vector class
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
                    assert(length(varargin) == 2, ...
                           'SReachTools:invalidArgs', ...
                           ['Gaussian random vector needs the mean vector ', ...
                           'and covariance matrix']);
                       
                    % Check if mean is a nonempty numeric column vector
                    validateattributes(varargin{1}, {'numeric'}, {'column'});
                    % Check if covariance matrix is a nonempty square matrix
                    validateattributes(varargin{2}, {'numeric'}, ...
                                       {'nonempty','square'});
                                   
                    % set parameters
                    obj.parameters.mean = varargin{1};
                    obj.parameters.covariance = varargin{2};
                    
                    % Check if the mean and covariance are of correct dimensions
                    assert(size(obj.parameters.mean,1) == ...
                           size(obj.parameters.covariance,1), ...
                           'SReachTools:invalidArgs', ...
                           ['Mean and covariance matrix have different ', ...
                            'dimensions']);
                        
                    % Update the dimension
                    obj.dim = size(obj.parameters.mean, 1);
                    
                    % Define an anonymous function using MATLAB's built-in pdf 
                    % Transpose the mean for mvnpdf
                    obj.pdf = @(x) mvnpdf(x, ...
                                          obj.parameters.mean', ...
                                          obj.parameters.covariance);
                otherwise
                    error('SReachTools:internal', ...
                          'Unsupported random vector type');
            end
        end
        
        function disp(obj)
        % SReachTools/RandomVector/disp  Override of MATLAB internal display
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
            
            disp(sprintf('%s-dimensional %s random vector', ...
                         num2str(obj.dim), ...
                         obj.type));
        end
    end
end
