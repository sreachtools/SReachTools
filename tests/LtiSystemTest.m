classdef LtiSystemTest < matlab.unittest.TestCase

    methods (Test)
        
        function testIncorrectEmptyFunctionCall(testCase)
            testCase.verifyError(@() LtiSystem(), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectNonSquareStateMatrix(testCase)
            testCase.verifyError(@() LtiSystem('StateMatrix',[1,1]), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectNoStateMatrixInput(testCase)
            testCase.verifyError(@() LtiSystem('InputMatrix', eye(2)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectInputMatrixStringOnly(testCase)
            testCase.verifyError(@() LtiSystem('InputMatrix'), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectInputMatrixBadString(testCase)
            testCase.verifyError(@() LtiSystem('InputMatrixGoneBad', eye(2)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectInputMatrixWrongRows(testCase)
            T = 0.5;
            testCase.verifyError(@() LtiSystem('StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2], ...
                'InputSpace', Polyhedron('lb', -1, 'ub', 1)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectInputMatrixWrongColumns(testCase)
            T = 0.5;
            testCase.verifyError(@() LtiSystem('StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2;T], ...
                'InputSpace', Polyhedron('lb', [-1;-1], 'ub', [1;1])), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectEmptyInputPolyhedronOneDimInputMatrix(testCase)
            T = 0.5;
            testCase.verifyError(@() LtiSystem('StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2;T]), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectInputPolyhedronOnly(testCase)
            T = 0.5;
            testCase.verifyError(@() LtiSystem('StateMatrix', [1, T; 0, 1], ...
                'Input', Polyhedron('lb', -1, 'ub', 1)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectNoStateMatrixDisturbance(testCase)
            testCase.verifyError(@() LtiSystem('DisturbanceMatrix', eye(2)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectDisturbanceMatrixStringOnly(testCase)
            testCase.verifyError(@() LtiSystem('DisturbanceMatrix'), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectDisturbanceMatrixBadString(testCase)
            testCase.verifyError(@() LtiSystem('DisturbanceMatrixGoneBad', eye(2)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectDisturbanceMatrixWrongRows(testCase)
            T = 0.5;
            testCase.verifyError(@() LtiSystem('StateMatrix', [1, T; 0, 1], ...
                'DisturbanceMatrix', [T^2], ...
                'DisturbanceSpace', Polyhedron('lb', -1, 'ub', 1)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectDisturbanceMatrixWrongColumns(testCase)
            T = 0.5;
            testCase.verifyError(@() LtiSystem('StateMatrix', [1, T; 0, 1], ...
                'DisturbanceMatrix', [T^2;T], ...
                'DisturbanceSpace', Polyhedron('lb', [-1;-1], 'ub', [1;1])), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectEmptyDisturbancePolyhedronOneDimDisturbanceMatrix(testCase)
            T = 0.5;
            testCase.verifyError(@() LtiSystem('StateMatrix', [1, T; 0, 1], ...
                'DisturbanceMatrix', [T^2;T]), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectStochasticDisturbanceBadDim(testCase)
            T = 0.5;
            mean_disturbance = zeros(5,1);
            covariance_disturbance = eye(5);
            GaussianDisturbance = RandomVector('Gaussian', ...
                mean_disturbance, ...
                covariance_disturbance);
            testCase.verifyError(@() LtiSystem('StateMatrix', [1, T; 0, 1], ...
                'DisturbanceMatrix', ones(2,4), ...
                'Disturbance', GaussianDisturbance), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectDisturbancePolyhedronOnly(testCase)
            T = 0.5;
            testCase.verifyError(@() LtiSystem('StateMatrix', [1, T; 0, 1], ...
                'Disturbance', Polyhedron('lb', -1, 'ub', 1)), 'SReachTools:invalidArgs');
        end
        
        function testCorrectInputPolyhedronOnly(testCase)
            T = 0.5;
            testCase.verifyClass(LtiSystem('StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2;T], ...
                'InputSpace', Polyhedron('lb', -1, 'ub', 1)), 'LtiSystem');
        end
        
        function testCorrectDisturbancePolyhedronOnly(testCase)
            T = 0.5;
            testCase.verifyClass(LtiSystem('StateMatrix', [1, T; 0, 1], ...
                'DisturbanceMatrix', [T^2;T], ...
                'Disturbance', Polyhedron('lb', -1, 'ub', 1)), 'LtiSystem');
        end

        function testCorrectArbitrarySystem(testCase)
            T = 0.5;
            testCase.verifyClass(LtiSystem('StateMatrix', zeros(2,2), ...
                'InputMatrix', ones(2,4), ...
                'InputSpace', Polyhedron('lb', -ones(4,1), 'ub', ones(4,1)), ...
                'DisturbanceMatrix', ones(2,6), ...
                'Disturbance', Polyhedron('lb', -ones(6,1), 'ub', ones(6,1))), 'LtiSystem');
        end

        function testCorrectDoubleIntegrator(testCase)
            T = 0.5;
            testCase.verifyClass(LtiSystem('StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2/2;T], ...
                'InputSpace', Polyhedron('lb', -1, 'ub', 1), ...
                'DisturbanceMatrix', [T^2/2;T], ...
                'Disturbance', Polyhedron('lb', -1, 'ub', 1)), 'LtiSystem');
        end
        
        function testCorrectDoubleIntegratorGaussian(testCase)
            T = 0.5;
            mean_disturbance = 0;
            covariance_disturbance = 1;
            GaussianDisturbance = RandomVector('Gaussian', ...
                mean_disturbance, ...
                covariance_disturbance);
            testCase.verifyClass(LtiSystem('StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2/2;T], ...
                'InputSpace', Polyhedron('lb', -1, 'ub', 1), ...
                'DisturbanceMatrix', [T^2/2;T], ...
                'Disturbance', GaussianDisturbance), 'LtiSystem');
        end

        function testCorrectDoubleIntegratorGaussianNoInput(testCase)
            T = 0.5;
            mean_disturbance = 0;
            covariance_disturbance = 1;
            GaussianDisturbance = RandomVector('Gaussian', ...
                mean_disturbance, ...
                covariance_disturbance);
            testCase.verifyClass(LtiSystem('StateMatrix', [1, T; 0, 1], ...
                'DisturbanceMatrix', [T^2/2;T], ...
                'Disturbance', GaussianDisturbance), 'LtiSystem');
        end
    end
end