classdef RandomVectorTests < matlab.unittest.TestCase
% SReachTools/TubeTests: Unit tests for bounded disturbances
% ===========================================================================
%
% Unit testing for RandomVector
%
% ===========================================================================
%
% This function is part of the Stochastic Optimal Control Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
%
%

    methods (Test)
        function inputHandlingTest(test_case)
            % Gaussian: Too many arguments => Varargin means invalidArgs
            test_case.verifyError(@(x) RandomVector('Gaussian',zeros(2,1),...
                eye(2), eye(2)), 'SReachTools:invalidArgs');
            % Gaussian: Not enough arguments
            test_case.verifyError(@(x) RandomVector('Gaussian',zeros(2,1)),...
                'SReachTools:invalidArgs');
            % Gaussian: Input parsing
            test_case.verifyError(@(x) RandomVector('Gaussian',zeros(2,1),...
                ones(2,3)), 'MATLAB:expectedSquare');
            test_case.verifyError(@(x) RandomVector('Gaussian',zeros(2,1),...
                ones(4,3)), 'MATLAB:expectedSquare');
            test_case.verifyError(@(x) RandomVector('Gaussian',zeros(1,3),...
                eye(3)), 'MATLAB:expectedColumn');
            % Gaussian: Dimension mismatch
            test_case.verifyError(@(x) RandomVector('Gaussian',zeros(2,1),...
                eye(3)), 'SReachTools:invalidArgs');
            % Gaussian: Non-symmetric matrix
            test_case.verifyWarning(@(x) RandomVector('Gaussian',zeros(2,1),...
                eye(2)+[0,3e-16;-3e-16,0]), 'SReachTools:runtime');
            RandomVector('Gaussian',zeros(2,1), eye(2)+[0,1e-18;-1e-18,0]);
            % Gaussian: Non-positive semi-definite matrix
            test_case.verifyError(@(x) RandomVector('Gaussian',zeros(2,1),...
                [-1,0;0,0]), 'SReachTools:invalidArgs');
            % Invalid string
            test_case.verifyError(@(x) RandomVector('Exp',zeros(2,1)),...
                'SReachTools:internal');
            % Simply define a well-defined Gaussian
            RandomVector('Gaussian',zeros(2,1), eye(2));
        end
        function multiplicationTest(test_case)        
            % Define a well-defined Gaussian
            r = RandomVector('Gaussian',zeros(2,1), eye(2));
            eye(2) * r;
            r * eye(2);
            [1, 0 ] * r;
            test_case.verifyWarning(@(x) [1, 0;1, 0;1, 0] * r,...
                'SReachTools:runtime');
            % Invalid dimension
            test_case.verifyError(@(x) r * eye(3), 'SReachTools:invalidArgs');            
            % Invalid input
            test_case.verifyError(@(x) r * 'ch', 'SReachTools:invalidArgs');            
        end
    end
end
