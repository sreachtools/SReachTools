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
%      https://github.com/abyvinod/Socbox/blob/master/LICENSE
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
        %      https://github.com/abyvinod/Socbox/blob/master/LICENSE
        %
        %

            test_case.verifyError(...
                @() StochasticDisturbance('Gaussian', zeros(5, 1)), ...
                'Socbox:invalidArgs', ...
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
        %      https://github.com/abyvinod/Socbox/blob/master/LICENSE
        %
        %
        
            test_case.verifyError(...
                @() StochasticDisturbance('Gaussian', ...
                    zeros(5, 1), eye(4)), ...
                'Socbox:invalidArgs', ...
                'Did not receive invalidArgs error');
        end

        function test1dGaussianPdfMatch(test_case)
        % Socbox/StochasticDisturbanceTests/test1dGaussianPdfMatch: Unit test 
        % for mistmatch in 1d Gaussian pdf
        % =====================================================================
        %
        % Unit test for mismatch in the StochasticDisturbance and known Gaussian
        % pdf for 1d Gaussian
        %
        % =====================================================================
        %
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://github.com/abyvinod/Socbox/blob/master/LICENSE
        %
        %

            disturbance = StochasticDisturbance('Gaussian', 1, 4);

            test_case.assertTrue(all(abs(normpdf([-10:0.1:10]', 1, 2) - ...
                disturbance.pdf([-10:0.1:10]')) < 1e-8 ), ...
                'Mismatch in the disturbance and known Gaussian pdfs');
        end

        function testNdGaussianPdfMatch(test_case)
        % Socbox/StochasticDisturbanceTests/test1dGaussianPdfMatch: Unit test 
        % for mistmatch in nd Gaussian pdf
        % =====================================================================
        %
        % Unit test for mismatch in the StochasticDisturbance and known Gaussian
        % pdf for nd Gaussian
        %
        % =====================================================================
        %
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://github.com/abyvinod/Socbox/blob/master/LICENSE
        %
        %
        
            dist_mean = zeros(5,1);
            dist_cov  = eye(5);
            sample_points = randn(100,5);

            disturbance = StochasticDisturbance('Gaussian', ...
                dist_mean, dist_cov);

            test_case.assertTrue(all(abs(...
                mvnpdf(sample_points, dist_mean', dist_cov) - ...
                disturbance.pdf(sample_points)) < 1e-8 ), ...
                'Mismatch in the disturbance and known Gaussian pdfs');

            test_case.assertTrue(disturbance.dimension == size(dist_cov, 2), ...
                ['Dimension mismatch between the disturbance object and ', ...
                 'the covariance matrix']);
        end
    end

end
