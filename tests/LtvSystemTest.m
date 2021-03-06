classdef LtvSystemTest < matlab.unittest.TestCase

    methods (Test)
        
        %% LtiSystem test copied with Lti replaced with Ltv
        function testIncorrectEmptyFunctionCall(testCase)
            testCase.verifyError(@() LtvSystem(), ...
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectNonSquareStateMatrix(testCase)
            testCase.verifyError(@() LtvSystem('StateMatrix',[1,1]), ...
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectNoStateMatrixInput(testCase)
            testCase.verifyError(@() LtvSystem('InputMatrix', eye(2)), ...
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectInputMatrixStringOnly(testCase)
            testCase.verifyError(@() LtvSystem('InputMatrix'), ...
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectInputMatrixBadString(testCase)
            testCase.verifyError(@() LtvSystem('InputMatrixGoneBad', ...
                eye(2)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectInputMatrixWrongRows(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2], ...
                'InputSpace', Polyhedron('lb', -1, 'ub', 1)), ...
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectInputMatrixWrongColumns(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2;T], ...
                'InputSpace', Polyhedron('lb', [-1;-1], 'ub', [1;1])), ...
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectEmptyInputPolyhedronOneDimInputMatrix(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2;T]), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectInputPolyhedronOnly(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', [1, T; 0, 1], ...
                'Input', Polyhedron('lb', -1, 'ub', 1)), ...
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectNoStateMatrixDisturbance(testCase)
            testCase.verifyError(@() LtvSystem('DisturbanceMatrix', ...
                eye(2)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectDisturbanceMatrixStringOnly(testCase)
            testCase.verifyError(@() LtvSystem('DisturbanceMatrix'), ...
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectDisturbanceMatrixBadString(testCase)
            testCase.verifyError(@() LtvSystem('DisturbanceMatrixGoneBad', ...
                eye(2)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectDisturbanceMatrixWrongRows(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', [1, T; 0, 1], ...
                'DisturbanceMatrix', [T^2], ...
                'DisturbanceSpace', Polyhedron('lb', -1, 'ub', 1)), ...
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectDisturbanceMatrixWrongColumns(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', [1, T; 0, 1], ...
                'DisturbanceMatrix', [T^2;T], ...
                'DisturbanceSpace', Polyhedron('lb', [-1;-1], 'ub', [1;1])), ...
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectEmptyDisturbancePolyhedronOneDimDisturbanceMatrix(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', [1, T; 0, 1], ...
                'DisturbanceMatrix', [T^2;T]), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectRandomVectorBadDim(testCase)
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
                'Disturbance', Polyhedron('lb', -1, 'ub', 1)), ...
                'SReachTools:invalidArgs');
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
                'Disturbance', Polyhedron('lb', -ones(6,1), 'ub', ones(6,1))), ...
                'LtvSystem');
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
            testCase.verifyError(@() LtvSystem('StateMatrix', @(t) [1,1]), ...
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectLtvNoStateMatrixInput(testCase)
            testCase.verifyError(@() LtvSystem('InputMatrix', @(t) eye(2)), ...
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectLtvInputMatrixBadString(testCase)
            testCase.verifyError(@() LtvSystem('InputMatrixGoneBad', ...
                @(t) eye(2)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectLtvInputMatrixWrongRows(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', @(t) [1, T; 0, 1], ...
                'InputMatrix', @(t) [T^2], ...
                'InputSpace', Polyhedron('lb', -1, 'ub', 1)), ...
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectLtvInputMatrixWrongColumns(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', @(t) [1, T; 0, 1], ...
                'InputMatrix', @(t) [T^2;T], ...
                'InputSpace', Polyhedron('lb', [-1;-1], 'ub', [1;1])), ...
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectLtvEmptyInputPolyhedronOneDimInputMatrix(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', ...
                @(t) [1, T; 0, 1], 'InputMatrix', @(t) [T^2;T]), ...
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectLtvInputPolyhedronOnly(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', ...
                @(t) [1, T; 0, 1], 'Input', Polyhedron('lb', -1, 'ub', 1)), ...
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectLtvNoStateMatrixDisturbance(testCase)
            testCase.verifyError(@() LtvSystem('DisturbanceMatrix', ...
                @(t) eye(2)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectLtvDisturbanceMatrixBadString(testCase)
            testCase.verifyError(@() LtvSystem('DisturbanceMatrixGoneBad', ...
                @(t) eye(2)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectLtvDisturbanceMatrixWrongRows(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', ...
                @(t) [1, T; 0, 1], ...
                'DisturbanceMatrix', @(t) [T^2], ...
                'DisturbanceSpace', Polyhedron('lb', -1, 'ub', 1)), ...
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectLtvDisturbanceMatrixWrongColumns(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', ...
                @(t) [1, T; 0, 1], ...
                'DisturbanceMatrix', @(t) [T^2;T], ...
                'DisturbanceSpace', Polyhedron('lb', [-1;-1], 'ub', [1;1])), ...
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectLtvEmptyDistPolyhedronOneDimDistMatrix(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', ...
                @(t) [1, T; 0, 1], ...
                'DisturbanceMatrix', @(t) [T^2;T]), ...
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectLtvRandomVectorBadDim(testCase)
            T = 0.5;
            mean_disturbance = zeros(5,1);
            covariance_disturbance = eye(5);
            GaussianDisturbance = RandomVector('Gaussian', ...
                mean_disturbance, ...
                covariance_disturbance);
            testCase.verifyError(@() LtvSystem('StateMatrix', ...
                @(t) [1, T; 0, 1], ...
                'DisturbanceMatrix', @(t) ones(2,4), ...
                'Disturbance', GaussianDisturbance), ...
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectLtvDisturbancePolyhedronOnly(testCase)
            T = 0.5;
            testCase.verifyError(@() LtvSystem('StateMatrix', ...
                @(t) [1, T; 0, 1], ...
                'Disturbance', Polyhedron('lb', -1, 'ub', 1)), ...
                'SReachTools:invalidArgs');
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
                'Disturbance', Polyhedron('lb', -ones(6,1), 'ub', ones(6,1))), ...
                'LtvSystem');
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
        
        function testInputHandlingLtvSystemMethods(testCase)
            T = 0.25;
            umax = 0.75;
            dmax = 1;
            mean_disturbance = 0;
            covariance_disturbance = 4;
            GaussianDisturbance = RandomVector('Gaussian', mean_disturbance, ...
                covariance_disturbance);
            sys = LtvSystem(...
                'StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2; T], ...
                'InputSpace', Polyhedron('lb', -umax, 'ub', umax), ...
                'DisturbanceMatrix', [T^2; T], ...
                'Disturbance', GaussianDisturbance);
            %% Tests getConcatMats,getConcatInputSpace
            %% Test for faulty time horizon
            initial_state = [1;1];
            % Gave a scalar negative time horizon
            time_horizon = -1;
%            testCase.verifyError(@() getHmatMeanCovForXSansInput(sys, ...
%                initial_state, time_horizon), 'SReachTools:invalidArgs');            
            testCase.verifyError(@() getConcatMats(sys, time_horizon),'SReachTools:invalidArgs');            
            testCase.verifyError(@() getConcatInputSpace(sys, time_horizon),'SReachTools:invalidArgs');
            % Gave a zero time horizon
            time_horizon = 0;
%            testCase.verifyError(@() getHmatMeanCovForXSansInput(sys, ...
%                initial_state, time_horizon), 'SReachTools:invalidArgs');            
            testCase.verifyError(@() getConcatMats(sys, time_horizon),'SReachTools:invalidArgs');            
            testCase.verifyError(@() getConcatInputSpace(sys, time_horizon),'SReachTools:invalidArgs');
            % Gave a vector of time horizon
            time_horizon = 0:10;
%            testCase.verifyError(@() getHmatMeanCovForXSansInput(sys, ...
%                [1;1], time_horizon), 'SReachTools:invalidArgs');
            testCase.verifyError(@() getConcatMats(sys, time_horizon),'SReachTools:invalidArgs');            
            testCase.verifyError(@() getConcatInputSpace(sys, time_horizon),'SReachTools:invalidArgs');
%            %% Test for faulty initial state
%            time_horizon = 10;
%            % Gave a row vector as an initial state
%            initial_state = [1 1];
%            testCase.verifyError(@() getHmatMeanCovForXSansInput(sys, ...
%                initial_state, time_horizon), 'SReachTools:invalidArgs');            
%            % Gave a matrix initial state
%            initial_state = eye(2);
%            testCase.verifyError(@() getHmatMeanCovForXSansInput(sys, ...
%                initial_state, time_horizon), 'SReachTools:invalidArgs');            
%            % Gave a random vector initial state with incorrect dimensions
%            initial_state = RandomVector('Gaussian',[1;1;1],eye(3));
%            testCase.verifyError(@() getHmatMeanCovForXSansInput(sys, ...
%                initial_state, time_horizon), 'SReachTools:invalidArgs');            
%            %% Test for non-stochastic system
%            sys = LtvSystem(...
%                'StateMatrix', [1, T; 0, 1], ...
%                'InputMatrix', [T^2; T], ...
%                'InputSpace', Polyhedron('lb', -umax, 'ub', umax), ...
%                'DisturbanceMatrix', [T^2; T], ...
%                'Disturbance', Polyhedron('lb', -dmax, 'ub', dmax));
%            % with non-stochastic initial state
%            initial_state = [1;1];
%            testCase.verifyError(@() getHmatMeanCovForXSansInput(sys, ...
%                initial_state, time_horizon), 'SReachTools:invalidArgs');            
%            % with random initial state
%            initial_state = RandomVector('Gaussian',[1;1;1],eye(3));
%            testCase.verifyError(@() getHmatMeanCovForXSansInput(sys, ...
%                initial_state, time_horizon), 'SReachTools:invalidArgs');            
            %% Test for uncontrolled system - getConcatInputSpace would fail
            time_horizon = 10;
            sys = LtvSystem(...
                'StateMatrix', [1, T; 0, 1], ...
                'DisturbanceMatrix', [T^2; T], ...
                'Disturbance', GaussianDisturbance);
            testCase.verifyError(@() getConcatInputSpace(sys, time_horizon),'SReachTools:invalidArgs');            
        end
        
%        function testgetHmatMeanCovForXSansInput(testCase)
%            time_horizon = 10;
%            T = 0.25;
%            umax = 0.75;
%            mean_disturbance = 0;
%            covariance_disturbance = 4;
%            GaussianDisturbance = RandomVector('Gaussian',mean_disturbance, ...
%                covariance_disturbance);                                                     
%            sys = LtvSystem(...
%                'StateMatrix', [1, T; 0, 1], ...
%                'InputMatrix', [T^2; T], ...
%                'InputSpace', Polyhedron('lb', -umax, 'ub', umax), ...
%                'DisturbanceMatrix', [T^2; T], ...
%                'Disturbance', GaussianDisturbance);
%            % Load Abar_saved, H_saved, G_matrix_saved expected
%            load('./data/getConcatMatsData.mat');
%            
%            %% Test for a non-stochastic initial state
%            initial_state = [2;0];
%            [H, mean_X_sans_input, cov_X_sans_input, Z, G] = ...
%               getHmatMeanCovForXSansInput(sys, ...
%                                           initial_state, ...
%                                           time_horizon);
%            testCase.verifyLessThanOrEqual(sum(sum(abs(Z - Abar_saved))),1e-8);
%            testCase.verifyLessThanOrEqual(sum(sum(abs(H - H_matrix_saved))),1e-8);
%            testCase.verifyLessThanOrEqual(sum(sum(abs(G - G_matrix_saved))),1e-8);
%            testCase.verifyLessThanOrEqual(sum(sum(abs(mean_X_sans_input - repmat(initial_state,time_horizon,1)))),1e-8);
%            testCase.verifyLessThanOrEqual(sum(sum(abs(...
%                cov_X_sans_input - covariance_disturbance * G_matrix_saved * G_matrix_saved'))), ...
%                1e-8);
%            
%            %% Test for a stochastic initial state
%            initial_state = RandomVector('Gaussian',[2;0],eye(2));
%            [H, mean_X_sans_input, cov_X_sans_input, Z, G] = ...
%               getHmatMeanCovForXSansInput(sys, ...
%                                           initial_state, ...
%                                           time_horizon);
%            testCase.verifyLessThanOrEqual(sum(sum(abs(Z - Abar_saved))),1e-8);
%            testCase.verifyLessThanOrEqual(sum(sum(abs(H - H_matrix_saved))),1e-8);
%            testCase.verifyLessThanOrEqual(sum(sum(abs(G - G_matrix_saved))),1e-8);
%            testCase.verifyLessThanOrEqual(sum(sum(abs(mean_X_sans_input - repmat([2;0],time_horizon,1)))),1e-8);
%            testCase.verifyLessThanOrEqual(sum(sum(abs(...
%                cov_X_sans_input - Abar_saved * Abar_saved'- covariance_disturbance * G_matrix_saved * G_matrix_saved'))), ...
%                1e-8);
%        end
        
        function testgetConcatMats(testCase)
            time_horizon = 10;
            T = 0.25;
            umax = 0.75;
            dmax = 1;
            mean_disturbance = 0;
            covariance_disturbance = 4;
            GaussianDisturbance = RandomVector('Gaussian',mean_disturbance, ...
                covariance_disturbance);                                                     
            % Load Abar_saved, H_saved, G_matrix_saved expected
            load('./data/getConcatMatsData.mat');
            
            sys = LtvSystem(...
                'StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2; T], ...
                'InputSpace', Polyhedron('lb', -umax, 'ub', umax), ...
                'DisturbanceMatrix', [T^2; T], ...
                'Disturbance', GaussianDisturbance);
            [Z,H,G] = getConcatMats(sys,time_horizon);
            testCase.verifyLessThanOrEqual(sum(sum(abs(Z - Abar_saved))),1e-8);
            testCase.verifyLessThanOrEqual(sum(sum(abs(H - H_matrix_saved))),1e-8);
            testCase.verifyLessThanOrEqual(sum(sum(abs(G - G_matrix_saved))),1e-8);
            % Gave a non-stochastic LtvSystem
            sys = LtvSystem(...
                'StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2; T], ...
                'InputSpace', Polyhedron('lb', -umax, 'ub', umax), ...
                'DisturbanceMatrix', [T^2; T], ...
                'Disturbance', Polyhedron('lb', -dmax, 'ub', dmax));
            [Z,H,G] = getConcatMats(sys,time_horizon);
            testCase.verifyLessThanOrEqual(sum(sum(abs(Z - Abar_saved))),1e-8);
            testCase.verifyLessThanOrEqual(sum(sum(abs(H - H_matrix_saved))),1e-8);
            testCase.verifyLessThanOrEqual(sum(sum(abs(G - G_matrix_saved))),1e-8);
            % Disturbance-free
            sys = LtvSystem(...
                'StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2; T], ...
                'InputSpace', Polyhedron('lb', -umax, 'ub', umax));
            [Z,H,G] = getConcatMats(sys,time_horizon);
            testCase.verifyLessThanOrEqual(sum(sum(abs(Z - Abar_saved))),1e-8);
            testCase.verifyLessThanOrEqual(sum(sum(abs(H - H_matrix_saved))),1e-8);
            testCase.verifyEqual(size(G),[2*time_horizon 0]);            
            % Control-free
            sys = LtvSystem(...
                'StateMatrix', [1, T; 0, 1], ...
                'DisturbanceMatrix', [T^2; T], ...
                'Disturbance', GaussianDisturbance);
            [Z,H,G] = getConcatMats(sys,time_horizon);
            testCase.verifyLessThanOrEqual(sum(sum(abs(Z - Abar_saved))),1e-8);
            testCase.verifyLessThanOrEqual(sum(sum(abs(G - G_matrix_saved))),1e-8);
            testCase.verifyEqual(size(H),[2*time_horizon 0]);            
        end       
        
        function testgetConcatInputSpace(testCase)
            time_horizon = 10;
            umax = 1;
            sys = LtvSystem('StateMatrix', eye(2), ...
                            'InputMatrix', ones(2,1), ...
                            'InputSpace', Polyhedron('lb', -umax, 'ub', umax));
            [concat_input_space_A, concat_input_space_b] = ...
                getConcatInputSpace(sys,time_horizon);
            obtained_polyhedron = Polyhedron('H',[concat_input_space_A, ...
                            concat_input_space_b]);
            expected_polyhedron = Polyhedron( ...
                'lb', -umax * ones(time_horizon,1), ...
                'ub',  umax * ones(time_horizon,1));
            testCase.verifyTrue(obtained_polyhedron == expected_polyhedron);
        end
        
        function testDisp(testCase)
        % Test display
        % We use evalc to skip stdout printing, but it will throw errors
        % We will use getChainOfIntegLtiSystem as the system generator
        % 
        % THIS MEANS WE CAN NOT TEST FOR WARNINGS
        %
        
            % Known turning rate sequence
            sampling_time = 0.1;                        % Sampling time
            time_horizon = 50;                          % Max number of time 
                                                        % steps in simulation
            init_heading = pi/10;                       % Initial heading 
            % Create a constant turning rate sequence
            omega = pi/time_horizon/sampling_time;      
            turning_rate = omega*ones(time_horizon,1);   
            % Input space definition
            umax = 6;
            input_space = Polyhedron('lb',0,'ub',umax);
            % Disturbance matrix and random vector definition
            dist_matrix = eye(2);
            eta_dist = RandomVector('Gaussian',zeros(2,1), 0.001 * eye(2));

            [sys, ~] = getDubinsCarLtv('add-dist', turning_rate, ...
              init_heading, sampling_time, input_space, dist_matrix, eta_dist);
            evalc('disp(sys)');
            evalc('disp(sys,''verbose'',true)');
            sysNoInpNoDist = LtvSystem('StateMatrix',sys.state_mat);
            evalc('disp(sysNoInpNoDist)');
            evalc('disp(sysNoInpNoDist,''verbose'',true)');
            sysNoDist = LtvSystem('StateMatrix',sys.state_mat, ...
                'InputMatrix',sys.input_mat, 'InputSpace', sys.input_space);
            evalc('disp(sysNoDist)');
            evalc('disp(sysNoDist,''verbose'',true)');
            sysNoInp = LtvSystem('StateMatrix',sys.state_mat, ...
                'DisturbanceMatrix',sys.dist_mat, 'Disturbance', sys.dist);
            evalc('disp(sysNoInp)');
            evalc('disp(sysNoInp,''verbose'',true)');
        end
    end
end
