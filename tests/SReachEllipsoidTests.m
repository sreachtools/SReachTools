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
            warning('on','SReachTools:runtime');
            test_case.verifyWarning(@(x) SReachEllipsoid(zeros(2,1),...
                eye(2)+[0,3e-16;-3e-16,0]), 'SReachTools:runtime');
            SReachEllipsoid(zeros(2,1), eye(2)+[0,1e-18;-1e-18,0]);
            % Gaussian: Non-positive semi-definite matrix
            test_case.verifyError(@(x) SReachEllipsoid(zeros(2,1),...
                [-1,0;0,0]), 'SReachTools:invalidArgs');
            % Simply define a well-defined ellipsoid -> circle
            SReachEllipsoid(zeros(2,1), eye(2));
        end
        
        function supportFunctionTest(test_case)
            % Support function of a unit circle must be one
            unit_cir = SReachEllipsoid(zeros(2,1), eye(2));
            theta = pi/4;
            v = [cos(theta), sin(theta)];
            test_case.verifyLessThanOrEqual(abs(unit_cir.support_fun(v)-1),...
                eps, 'All tangents of unit circle must return 1');
            theta = 3*pi/4;
            v = [cos(theta), sin(theta)];
            test_case.verifyLessThanOrEqual(abs(unit_cir.support_fun(v)-1),...
                eps, 'All tangents of unit circle must return 1');
            
            % Checking multiple support function vector candidates +
            % Cross-validate with CVX results
            rand_ellipse = SReachEllipsoid(ones(2,1), 10*eye(2)+[1,0;0,5]);
            theta = 3*pi/4;
            v = [cos(theta), sin(theta);
                 cos(1.2*pi/2+theta), sin(1.2*pi/2+theta)];
            support_fun_vals = rand_ellipse.support_fun(v);
            Q_inv_sqrt = chol(inv(rand_ellipse.shape_matrix));
            cvx_begin quiet
                variable y1(2,1);
                maximize (v(1,:)*y1)
                subject to
                    norm(Q_inv_sqrt * (y1 - rand_ellipse.center)) <= 1
            cvx_end
            cvx_begin quiet
                variable y2(2,1);
                maximize (v(2,:)*y2)
                subject to
                    norm(Q_inv_sqrt * (y2 - rand_ellipse.center)) <= 1
            cvx_end
                        
            test_case.verifyLessThanOrEqual(...
                abs(support_fun_vals- diag(v*[y1,y2])),...
                1e-8, 'CVX and support function definition do not agree!');
        end
        
        function multiplyTest(test_case)
            % Support function of a unit circle must be one
            unit_cir = SReachEllipsoid(ones(2,1), eye(2));
            
            F = [1,0;
                 0,3];
            
            ellipse_from_cir = F * unit_cir;
            
            % Square-root of eigenvalues is the semi-axis length
            test_case.verifyTrue(isequal(eig(ellipse_from_cir.shape_matrix),...
                [1;9]),...
                'Incorrect axis lengths of the ellipse constructed by scaling');
            test_case.verifyTrue(isequal(ellipse_from_cir.center,[1;3]),...
                'Incorrect center of the ellipse constructed by scaling');
            
        end
    end
end
