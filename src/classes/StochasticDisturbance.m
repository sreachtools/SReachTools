classdef StochasticDisturbance < RandomVector
% SReachTools/StochasticDisturbance: Create a stochastic disturbance object
% (subclass of RandomVector)
% ==========================================================================
%
% Defines a stochastic disturbance with a standard probability density function
% (pdf)
%
% We support the following pdfs:
%     1. Gaussian
%
% Usage
% ------
%
% Define a Gaussian stochastic disturbance of mean 0 and standard deviation 2:
% GaussianRV = StochasticDisturbance('Gaussian', ...
%                                    0,
%                                    2^2);
%
% Define a Gaussian stochastic disturbance of mean [0;2] and covariance matrix 
% eye(2):
% GaussianRV = StochasticDisturbance('Gaussian', ...
%                                    [0;2],
%                                    eye(2));
%   
% ==========================================================================
%
% StochasticDisturbance Properties:
% ------------------------
%   type       - Stochastic disturbance type (string)
%   parameters - System parameters (struct)
%   dim        - Stochastic disturbance dimension (scalar)
%   pdf        - Probability density function (function handle)
%
% StochasticDisturbance Methods:
% ---------------------
%   StochasticDisturbance/StochasticDisturbance - Class constructor
% 
% Notes:
% ------
% * MATLAB DEPENDENCY: Uses MATLAB's Statistics and Machine Learning Toolbox.
%                      Needs mvnpdf
% * Currently only supports Gaussian stochastic disturbance
% * Requires the mean and the covariance matrices to be non-empty column
%   vector and a symmetric matrix respectively
% * The anonymous function used for the definition of obj.pdf transposes the
%   accepted column vector for using mvnpdf.
% * StochasticDisturbance.pdf takes in arguments of the form N_points x
%   stochastic_disturbance_dim
% 
% =========================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 
% 
    methods
        function obj = StochasticDisturbance(rv_type, varargin)
            switch(lower(rv_type))
                case 'gaussian'
                    if length(varargin) ~= 2
                        throwAsCaller(SrtInvalidArsError(['Gaussian disturbance needs the mean vector ', ...
                           'and covariance matrix']));
                    end
                    
                    args{1} = varargin{1};
                    args{2} = varargin{2};
                otherwise
                    exc = SrtInvalidArgsError('Unsupported disturbance type');
                    throw(exc);
            end
            % Call to superclass constructor must be unconditional
            obj = obj@RandomVector(rv_type, args{:});
        end

        function disp(obj)
        % SReachTools/StochasticDisturbance/disp  Override of MATLAB internal
        % display
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
            
            disp(sprintf('%s-dimensional %s stochastic disturbance', ...
                         num2str(obj.dim), ...
                         obj.type));
        end
    end
end
