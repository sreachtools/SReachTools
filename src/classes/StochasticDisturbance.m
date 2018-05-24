classdef StochasticDisturbance
% SReachTools/StochasticDisturbance: Create a stochastic disturbance object
% ==========================================================================
%
% Defines a stochastic disturbance with a standard probability density
% function (pdf)
%
% We support the following pdfs:
%     1. Gaussian
%
% Usage:
% ------
%
% % Define a Gaussian disturbance of mean 0 and standard deviation 2:
% GaussianDisturbance = STOCHASTICDISTURBANCE('Gaussian',...
%                                             0,
%                                             2^2);
%
% % Define a Gaussian disturbance of mean [0;2] and covariance matrix eye(2):
% GaussianDisturbance = STOCHASTICDISTURBANCE('Gaussian',...
%                                             [0;2],
%                                             eye(2));
%   
% ==========================================================================
%
% STOCHASTICDISTURBANCE Properties:
% ---------------------------------
%   type       - Disturbance type (string)
%   parameters - System parameters (struct)
%   dimension  - Disturbance dimension (scalar)
%   pdf        - Probability density function (function handle)
%
% STOCHASTICDISTURBANCE Methods:
% ------------------------------
%   StochasticDisturbance/StochasticDisturbance - Class constructor
% 
% Notes:
% ------
% * MATLAB DEPENDENCY: Uses MATLAB's Statistics and Machine Learning Toolbox.
%                      Needs mvnpdf
% * Currently only supports Gaussian Distrubances
% * Requires the mean and the covariance matrices to be non-empty column
%   vector and a symmetric matrix respectively
% * mvnpdf requires mean to be a row vector. Hence, while the class accepts a
%   column vector for the mean of a Gaussian disturbance, the anonymous
%   function definition used for obj.pdf transposes it.
% * StochasticDisturbance.pdf takes in arguments of the form 
%   N_points x disturbance_dimension
% 
% =========================================================================
% 
% This function is part of the Stochastic Optimal Control Toolbox.
% License for the use of this function is given in
%      https://github.com/abyvinod/SReachTools/blob/master/LICENSE
% 
% 

    properties
        % StochasticDisturbance/type
        % ==================================================================
        % 
        % String indicator of type of disturbance
        %
        % Acceptable types:
        %   Gaussian
        % 
        type

        % StochasticDisturbance/parameters
        % ==================================================================
        % 
        % Struct containing stochastic disturbance parameter information; 
        % will vary for different disturbance types. Currently only Gaussian 
        % disturbances are supported.
        %
        % Gaussian type:
        %   parameters.mean       - Mean vector (n x 1)
        %   parameters.covariance - Covariance matrix (n x n)
        % 
        parameters
        dimension
        pdf
    end
    methods
        function obj = StochasticDisturbance(disturbance_type, varargin)
        % SReachTools/StochasticDisturbance: Constructor for StochasticDisturbance
        % class
        % ====================================================================
        %
        % Stochastic Disturbance class constructor
        %
        % Usage
        % ------
        % % Define a Gaussian disturbance of mean 0 and standard deviation 2:
        % GaussianDisturbance = STOCHASTICDISTURBANCE('Gaussian',...
        %                                             0,
        %                                             2^2);
        %
        % % Define a Gaussian disturbance of mean [0;2] and covariance matrix 
        % % eye(2):
        % GaussianDisturbance = STOCHASTICDISTURBANCE('Gaussian',...
        %                                             [0;2],
        %                                             eye(2));
        %
        % =====================================================================
        %
        % Inputs:
        % -------
        %   disturbance_type - Character array of disturbance type; currently
        %                      acceptable types:
        %                          'Gaussian'
        %   varargin         - Arguments with the properties of the disturbance
        %                      For Gaussian:
        %                          mu
        %                          sigma
        %
        % Ouptus:
        % -------
        %   obj - Stochastic Disturbance object
        %
        % =====================================================================
        % 
        %   This function is part of the Stochastic Optimal Control Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/abyvinod/SReachTools/blob/master/LICENSE
        % 
        % 
            
            % Check if the disturbance type is a string 
            validateattributes(disturbance_type, {'char'}, {'nonempty'});
            % Update the stochastic disturbance type
            obj.type = disturbance_type; 

            switch(lower(obj.type))
                case 'gaussian'
                    assert(length(varargin) == 2,...
                           'SReachTools:invalidArgs',...
                           ['Gaussian disturbance needs the mean vector and ',...
                           'covariance matrix']);
                       
                    % Check if mean is a nonempty numeric column vector
                    validateattributes(varargin{1}, {'numeric'}, {'column'});
                    % Check if covariance matrix is a nonempty square matrix
                    validateattributes(varargin{2}, {'numeric'},...
                                       {'nonempty','square'});
                                   
                    % set parameters
                    obj.parameters.mean = varargin{1};
                    obj.parameters.covariance = varargin{2};
                    
                    % Check if the mean and covariance are of correct dimensions
                    assert(size(obj.parameters.mean,1) ==...
                           size(obj.parameters.covariance,1),...
                           'SReachTools:invalidArgs',...
                           ['Mean and covariance matrix have different ',...
                            'dimensions']);
                        
                    % Update the dimension
                    obj.dimension = size(obj.parameters.mean, 1);
                    
                    % Define an anonymous function using MATLAB's built-in pdf 
                    % Transpose the mean for mvnpdf
                    obj.pdf = @(x) mvnpdf(x,...
                                          obj.parameters.mean',...
                                          obj.parameters.covariance);
                otherwise
                    error('SReachTools:internal',...
                          'Unsupported disturbance type');
            end
        end
        
        function disp(obj)
        % SReachTools/StochasticDisturbance/disp  Override of MATLAB internal display
        % ====================================================================
        % 
        % Overriding of MATLAB built-in display function for the class
        %
        % Usage:
        % % omitted
        %
        % ====================================================================
        %
        % disp(obj)
        %
        % Inputs:  None
        % Outputs: None
        % 
        % ====================================================================
        % 
        %   This function is part of the Stochastic Optimal Control Toolbox.
        %   License for the use of this function is given in
        %        https://github.com/abyvinod/SReachTools/blob/master/LICENSE
        % 
        %
            
            disp(sprintf('%s-dimensional %s disturbance',...
                         num2str(obj.dimension),...
                         obj.type));
        end
    end
end
