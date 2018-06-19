classdef LtvSystemTest < matlab.unittest.TestCase

    methods (Test)
        
        %% LtiSystem test copied with Lti replaced with Ltv
        function testIncorrectEmptyFunctionCall(testCase)
            testCase.verifyError(@() LtvSystem(), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectNonSquareStateMatrix(testCase)
            testCase.verifyError(@() LtvSystem('StateMatrix',[1,1]), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectNoStateMatrixInput(testCase)
            testCase.verifyError(@() LtvSystem('InputMatrix', eye(2)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectInputMatrixStringOnly(testCase)
            testCase.verifyError(@() LtvSystem('InputMatrix'), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectInputMatrixBadString(testCase)
            testCase.verifyError(@() LtvSystem('InputMatrixGoneBad', eye(2)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectInputMatrixWrongRows(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2], ...
                'InputSpace', Polyhedron('lb', -1, 'ub', 1)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectInputMatrixWrongColumns(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2;T], ...
                'InputSpace', Polyhedron('lb', [-1;-1], 'ub', [1;1])), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectEmptyInputPolyhedronOneDimInputMatrix(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2;T]), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectInputPolyhedronOnly(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', [1, T; 0, 1], ...
                'Input', Polyhedron('lb', -1, 'ub', 1)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectNoStateMatrixDisturbance(testCase)
            testCase.verifyError(@() LtvSystem('DisturbanceMatrix', eye(2)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectDisturbanceMatrixStringOnly(testCase)
            testCase.verifyError(@() LtvSystem('DisturbanceMatrix'), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectDisturbanceMatrixBadString(testCase)
            testCase.verifyError(@() LtvSystem('DisturbanceMatrixGoneBad', eye(2)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectDisturbanceMatrixWrongRows(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', [1, T; 0, 1], ...
                'DisturbanceMatrix', [T^2], ...
                'DisturbanceSpace', Polyhedron('lb', -1, 'ub', 1)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectDisturbanceMatrixWrongColumns(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', [1, T; 0, 1], ...
                'DisturbanceMatrix', [T^2;T], ...
                'DisturbanceSpace', Polyhedron('lb', [-1;-1], 'ub', [1;1])), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectEmptyDisturbancePolyhedronOneDimDisturbanceMatrix(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', [1, T; 0, 1], ...
                'DisturbanceMatrix', [T^2;T]), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectStochasticDisturbanceBadDim(testCase)
            T = 0.5;
            mean_disturbance = zeros(5,1);
            covariance_disturbance = eye(5);
            GaussianDisturbance = RandomVector('Gaussian', ...
                mean_disturbance, ...
                covariance_disturbance);
            testCase.verifyError(@() LtvSystem('StateMatrix', [1, T; 0, 1], ...
                'DisturbanceMatrix', ones(2,4), ...
                'Disturbance', GaussianDisturbance), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectDisturbancePolyhedronOnly(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', [1, T; 0, 1], ...
                'Disturbance', Polyhedron('lb', -1, 'ub', 1)), 'SReachTools:invalidArgs');
        end
        
        function testCorrectInputPolyhedronOnly(testCase)
            T = 0.5;
            testCase.verifyClass(LtvSystem('StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2;T], ...
                'InputSpace', Polyhedron('lb', -1, 'ub', 1)), 'LtvSystem');
        end
        
        function testCorrectDisturbancePolyhedronOnly(testCase)
            T = 0.5;
            testCase.verifyClass(LtvSystem('StateMatrix', [1, T; 0, 1], ...
                'DisturbanceMatrix', [T^2;T], ...
                'Disturbance', Polyhedron('lb', -1, 'ub', 1)), 'LtvSystem');
        end

        function testCorrectArbitrarySystem(testCase)
            T = 0.5;
            testCase.verifyClass(LtvSystem('StateMatrix', zeros(2,2), ...
                'InputMatrix', ones(2,4), ...
                'InputSpace', Polyhedron('lb', -ones(4,1), 'ub', ones(4,1)), ...
                'DisturbanceMatrix', ones(2,6), ...
                'Disturbance', Polyhedron('lb', -ones(6,1), 'ub', ones(6,1))), 'LtvSystem');
        end

        function testCorrectDoubleIntegrator(testCase)
            T = 0.5;
            testCase.verifyClass(LtvSystem('StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2/2;T], ...
                'InputSpace', Polyhedron('lb', -1, 'ub', 1), ...
                'DisturbanceMatrix', [T^2/2;T], ...
                'Disturbance', Polyhedron('lb', -1, 'ub', 1)), 'LtvSystem');
        end
        
        function testCorrectDoubleIntegratorGaussian(testCase)
            T = 0.5;
            mean_disturbance = 0;
            covariance_disturbance = 1;
            GaussianDisturbance = RandomVector('Gaussian', ...
                mean_disturbance, ...
                covariance_disturbance);
            testCase.verifyClass(LtvSystem('StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2/2;T], ...
                'InputSpace', Polyhedron('lb', -1, 'ub', 1), ...
                'DisturbanceMatrix', [T^2/2;T], ...
                'Disturbance', GaussianDisturbance), 'LtvSystem');
        end

        function testCorrectDoubleIntegratorGaussianNoInput(testCase)
            T = 0.5;
            mean_disturbance = 0;
            covariance_disturbance = 1;
            GaussianDisturbance = RandomVector('Gaussian', ...
                mean_disturbance, ...
                covariance_disturbance);
            testCase.verifyClass(LtvSystem('StateMatrix', [1, T; 0, 1], ...
                'DisturbanceMatrix', [T^2/2;T], ...
                'Disturbance', GaussianDisturbance), 'LtvSystem');
        end
        
        %% LtvSystem novel tests
        function testIncorrectLtvNonSquareStateMatrix(testCase)
            testCase.verifyError(@() LtvSystem('StateMatrix', @(t) [1,1]), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectLtvNoStateMatrixInput(testCase)
            testCase.verifyError(@() LtvSystem('InputMatrix', @(t) eye(2)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectLtvInputMatrixBadString(testCase)
            testCase.verifyError(@() LtvSystem('InputMatrixGoneBad', @(t) eye(2)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectLtvInputMatrixWrongRows(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', @(t) [1, T; 0, 1], ...
                'InputMatrix', @(t) [T^2], ...
                'InputSpace', Polyhedron('lb', -1, 'ub', 1)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectLtvInputMatrixWrongColumns(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', @(t) [1, T; 0, 1], ...
                'InputMatrix', @(t) [T^2;T], ...
                'InputSpace', Polyhedron('lb', [-1;-1], 'ub', [1;1])), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectLtvEmptyInputPolyhedronOneDimInputMatrix(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', @(t) [1, T; 0, 1], ...
                'InputMatrix', @(t) [T^2;T]), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectLtvInputPolyhedronOnly(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', @(t) [1, T; 0, 1], ...
                'Input', Polyhedron('lb', -1, 'ub', 1)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectLtvNoStateMatrixDisturbance(testCase)
            testCase.verifyError(@() LtvSystem('DisturbanceMatrix', @(t) eye(2)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectLtvDisturbanceMatrixBadString(testCase)
            testCase.verifyError(@() LtvSystem('DisturbanceMatrixGoneBad', @(t) eye(2)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectLtvDisturbanceMatrixWrongRows(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', @(t) [1, T; 0, 1], ...
                'DisturbanceMatrix', @(t) [T^2], ...
                'DisturbanceSpace', Polyhedron('lb', -1, 'ub', 1)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectLtvDisturbanceMatrixWrongColumns(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', @(t) [1, T; 0, 1], ...
                'DisturbanceMatrix', @(t) [T^2;T], ...
                'DisturbanceSpace', Polyhedron('lb', [-1;-1], 'ub', [1;1])), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectLtvEmptyDistPolyhedronOneDimDistMatrix(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', @(t) [1, T; 0, 1], ...
                'DisturbanceMatrix', @(t) [T^2;T]), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectLtvStochasticDisturbanceBadDim(testCase)
            T = 0.5;
            mean_disturbance = zeros(5,1);
            covariance_disturbance = eye(5);
            GaussianDisturbance = RandomVector('Gaussian', ...
                mean_disturbance, ...
                covariance_disturbance);
            testCase.verifyError(@() LtvSystem('StateMatrix', @(t) [1, T; 0, 1], ...
                'DisturbanceMatrix', @(t) ones(2,4), ...
                'Disturbance', GaussianDisturbance), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectLtvDisturbancePolyhedronOnly(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', @(t) [1, T; 0, 1], ...
                'Disturbance', Polyhedron('lb', -1, 'ub', 1)), 'SReachTools:invalidArgs');
        end
        
        function testCorrectLtvInputPolyhedronOnly(testCase)
            T = 0.5;
            testCase.verifyClass(LtvSystem('StateMatrix', @(t) [1, T; 0, 1], ...
                'InputMatrix', @(t) [T^2;T], ...
                'InputSpace', Polyhedron('lb', -1, 'ub', 1)), 'LtvSystem');
        end
        
        function testCorrectLtvDisturbancePolyhedronOnly(testCase)
            T = 0.5;
            testCase.verifyClass(LtvSystem('StateMatrix', @(t) [1, T; 0, 1], ...
                'DisturbanceMatrix', @(t) [T^2;T], ...
                'Disturbance', Polyhedron('lb', -1, 'ub', 1)), 'LtvSystem');
        end

        function testCorrectLtvArbitrarySystem(testCase)
            T = 0.5;
            testCase.verifyClass(LtvSystem('StateMatrix', @(t) zeros(2,2), ...
                'InputMatrix', @(t) ones(2,4), ...
                'InputSpace', Polyhedron('lb', -ones(4,1), 'ub', ones(4,1)), ...
                'DisturbanceMatrix', @(t) ones(2,6), ...
                'Disturbance', Polyhedron('lb', -ones(6,1), 'ub', ones(6,1))), 'LtvSystem');
        end

        function testCorrectLtvDoubleIntegrator(testCase)
            T = 0.5;
            testCase.verifyClass(LtvSystem('StateMatrix', @(t) [1, T; 0, 1], ...
                'InputMatrix', @(t) [T^2/2;T], ...
                'InputSpace', Polyhedron('lb', -1, 'ub', 1), ...
                'DisturbanceMatrix', @(t) [T^2/2;T], ...
                'Disturbance', Polyhedron('lb', -1, 'ub', 1)), 'LtvSystem');
        end
        
        function testCorrectLtvDoubleIntegratorGaussian(testCase)
            T = 0.5;
            mean_disturbance = 0;
            covariance_disturbance = 1;
            GaussianDisturbance = RandomVector('Gaussian', ...
                mean_disturbance, ...
                covariance_disturbance);
            testCase.verifyClass(LtvSystem('StateMatrix', @(t) [1, T; 0, 1], ...
                'InputMatrix', @(t) [T^2/2;T], ...
                'InputSpace', Polyhedron('lb', -1, 'ub', 1), ...
                'DisturbanceMatrix', @(t) [T^2/2;T], ...
                'Disturbance', GaussianDisturbance), 'LtvSystem');
        end

        function testCorrectLtvDoubleIntegratorGaussianNoInput(testCase)
            T = 0.5;
            mean_disturbance = 0;
            covariance_disturbance = 1;
            GaussianDisturbance = RandomVector('Gaussian', ...
                mean_disturbance, ...
                covariance_disturbance);
            testCase.verifyClass(LtvSystem('StateMatrix', @(t) [1, T; 0, 1], ...
                'DisturbanceMatrix', @(t) [T^2/2;T], ...
                'Disturbance', GaussianDisturbance), 'LtvSystem');
        end
    end
end