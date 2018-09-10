classdef StochasticDisturbanceTests < matlab.unittest.TestCase
% Socbox/StochasticDisturbanceTests: Unit tests for bounded disturbances
% ===========================================================================
%
% Unit testing for bounded disturbances
%
% Usage:
% ------
% tests = StochasticDisturbanceTests()
% run(tests)
%
% ===========================================================================
%
% This function is part of the Stochastic Optimal Control Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/Socbox/blob/master/LICENSE
%
%

    methods (Test)
        function testIncompleteGaussianDisturbance(test_case)
        % Socbox/StochasticDisturbanceTests/testIncompleteGaussianDisturbance: 
        % Unit test for imcompletely specified Gaussian disturbance
        % =====================================================================
        %
        % Unit test for imcompletely specified Gaussian disturbance
        %
        % Should throw a 'Socbox:invalidArgs' error
        %
        % =====================================================================
        %
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://github.com/unm-hscl/Socbox/blob/master/LICENSE
        %
        %

            test_case.verifyError(...
                @() StochasticDisturbance('Gaussian', zeros(5, 1)), ...
                ?SrtInvalidArgsError, ...
                'Did not receive invalidArgs error');
        end

        function testGaussianDisturbanceDimensionMismatch(test_case)
        % Socbox/StochasticDisturbanceTests/testGaussianDisturbanceDimensionMismatch: 
        % Unit test for mistmatch in Gaussian disturbance
        % =====================================================================
        %
        % Unit test for mistmatch in Gaussian disturbance, i.e. mean and
        % covariance are not equally sized
        %
        % Should throw a 'Socbox:invalidArgs' error
        %
        % =====================================================================
        %
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://github.com/unm-hscl/Socbox/blob/master/LICENSE
        %
        %
        
            test_case.verifyError(...
                @() StochasticDisturbance('Gaussian', ...
                    zeros(5, 1), eye(4)), ...
                ?SrtInvalidArgsError, ...
                'Did not receive invalidArgs error');
        end
    end

end
