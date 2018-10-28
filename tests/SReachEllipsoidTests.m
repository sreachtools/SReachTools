classdef SReachEllipsoidTests < matlab.unittest.TestCase
% SReachTools/TubeTests: Unit tests for bounded disturbances
% ===========================================================================
%
% Unit testing for SReachEllipsoid
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
            % Gaussian: Too many arguments
            test_case.verifyError(@(x) SReachEllipsoid(zeros(2,1),...
                eye(2), eye(2)), 'MATLAB:TooManyInputs');
            % Gaussian: Not enough arguments
            test_case.verifyError(@(x) SReachEllipsoid(zeros(2,1)),...
                'SReachTools:invalidArgs');
            % Gaussian: Input parsing
            test_case.verifyError(@(x) SReachEllipsoid(zeros(2,1),...
                ones(2,3)), 'SReachTools:invalidArgs');
            test_case.verifyError(@(x) SReachEllipsoid(zeros(2,1),...
                ones(4,3)), 'SReachTools:invalidArgs');
            test_case.verifyError(@(x) SReachEllipsoid(zeros(1,3),...
                eye(3)), 'SReachTools:invalidArgs');
            % Gaussian: Dimension mismatch
            test_case.verifyError(@(x) SReachEllipsoid(zeros(2,1),...
                eye(3)), 'SReachTools:invalidArgs');
            % Gaussian: Non-symmetric matrix
            test_case.verifyWarning(@(x) SReachEllipsoid(zeros(2,1),...
                eye(2)+[0,3e-16;-3e-16,0]), 'SReachTools:runtime');
            SReachEllipsoid(zeros(2,1), eye(2)+[0,1e-18;-1e-18,0]);
            % Gaussian: Non-positive semi-definite matrix
            test_case.verifyError(@(x) SReachEllipsoid(zeros(2,1),...
                [-1,0;0,0]), 'SReachTools:invalidArgs');
            % Simply define a well-defined ellipsoid -> circle
            SReachEllipsoid(zeros(2,1), eye(2));
        end
    end
end
