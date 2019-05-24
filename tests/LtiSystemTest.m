classdef LtiSystemTest < matlab.unittest.TestCase
% LtiSystemTest
% ============================================================================
% 
% Unit tests for LtiSystem class
% 
% ============================================================================
% 
% Tests:
% ------
%   testIncorrectEmptyFunctionCall
%   testIncorrectNonSquareStateMatrix
%   testIncorrectNoStateMatrixInput
%   testIncorrectInputMatrixStringOnly
%   testIncorrectInputMatrixBadString
%   testIncorrectInputMatrixWrongRows
%   testIncorrectInputMatrixWrongColumns
%   testIncorrectEmptyInputPolyhedronOneDimInputMatrix
%   testIncorrectInputPolyhedronOnly
%   testIncorrectNoStateMatrixDisturbance
%   testIncorrectDisturbanceMatrixStringOnly
%   testIncorrectDisturbanceMatrixBadString
%   testIncorrectDisturbanceMatrixWrongRows
%   testIncorrectDisturbanceMatrixWrongColumns
%   testIncorrectEmptyDisturbancePolyhedronOneDimDisturbanceMatrix
%   testIncorrectRandomVectorBadDim
%   testIncorrectDisturbancePolyhedronOnly
%   testCorrectInputPolyhedronOnly
%   testCorrectDisturbancePolyhedronOnly
%   testCorrectArbitrarySystem
%   testCorrectDoubleIntegrator
%   testCorrectDoubleIntegratorGaussian
%   testCorrectDoubleIntegratorGaussianNoInput
%   testgetConcatMats
%   testgetConcatInputSpace
%   testDisp
%   
% ============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://sreachtools.github.io/license/
% 
% 

    methods (Test)
        
        function testIncorrectEmptyFunctionCall(testCase)
        % Test for failure when calling LtiSystem without inputs
        % 
        % Expected Error: SReachTools:invalidArgs
        % 

            testCase.verifyError(@() LtiSystem(), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectNonSquareStateMatrix(testCase)
        % Test for failure when calling LtiSystem with non-square state matrix
        % 
        % Expected Error: SReachTools:invalidArgs
        % 

            testCase.verifyError(@() LtiSystem('StateMatrix',[1,1]), ...
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectNoStateMatrixInput(testCase)
        % Test for failure when calling LtiSystem with no state matrix
        % 
        % Expected Error: SReachTools:invalidArgs
        % 

            testCase.verifyError(@() LtiSystem('InputMatrix', eye(2)), ....
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectInputMatrixStringOnly(testCase)
        % Test for failure when calling LtiSystem with only 'InputMatrix' string
        % 
        % Expected Error: SReachTools:invalidArgs
        % 
        
            testCase.verifyError(@() LtiSystem('InputMatrix'), ...
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectInputMatrixBadString(testCase)
        % Test for failure when calling LtiSystem with invalid property string
        % 
        % Expected Error: SReachTools:invalidArgs
        % 
        
            testCase.verifyError(@() LtiSystem('InputMatrixGoneBad', ...
                eye(2)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectInputMatrixWrongRows(testCase)
        % Test for failure when calling LtiSystem with mismatched State and 
        % input matrix dimension
        % 
        % Expected Error: SReachTools:invalidArgs
        % 
        
            T = 0.5;
            testCase.verifyError(@() LtiSystem('StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2], ...
                'InputSpace', Polyhedron('lb', -1, 'ub', 1)), ...
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectInputMatrixWrongColumns(testCase)
        % Test for failure when calling LtiSystem with mistmatch in input matrix
        % dimension and input space dimension
        % 
        % Expected Error: SReachTools:invalidArgs
        % 
        
            T = 0.5;
            testCase.verifyError(@() LtiSystem('StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2;T], ...
                'InputSpace', Polyhedron('lb', [-1;-1], 'ub', [1;1])), ...
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectEmptyInputPolyhedronOneDimInputMatrix(testCase)
        % Test for failure when calling LtiSystem with non-zero input matrix
        % and no defined input space
        % 
        % Expected Error: SReachTools:invalidArgs
        % 
        
            T = 0.5;
            testCase.verifyError(@() LtiSystem('StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2;T]), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectInputPolyhedronOnly(testCase)
        % Test for failure when calling LtiSystem with invalid input space 
        % string
        % 
        % Expected Error: SReachTools:invalidArgs
        % 
        % TODO: This test is basically the same as 
        %       testIncorrectInputMatrixBadString because they both have invalid
        %       input strings; should combine
        
            T = 0.5;
            testCase.verifyError(@() LtiSystem('StateMatrix', [1, T; 0, 1], ...
                'Input', Polyhedron('lb', -1, 'ub', 1)), ...
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectNoStateMatrixDisturbance(testCase)
        % Test for failure when calling LtiSystem with disturbance matrix but
        % no state matrix
        % 
        % Expected Error: SReachTools:invalidArgs
        % 
        % TODO: This is also the same as the test in which we specify empty /
        %       input matrix but no state matrix; combine
        % 

            testCase.verifyError(@() LtiSystem('DisturbanceMatrix', ...
                eye(2)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectDisturbanceMatrixStringOnly(testCase)
        % Test for failure when calling LtiSystem with disturbance matrix string
        % but no matrix
        % 
        % Expected Error: SReachTools:invalidArgs
        % 

            testCase.verifyError(@() LtiSystem('DisturbanceMatrix'), ...
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectDisturbanceMatrixBadString(testCase)
        % Test for failure when calling LtiSystem with bas disturbance matrix 
        % string
        % 
        % Expected Error: SReachTools:invalidArgs
        % 
        % TODO: Again, repeated test for the bad input matrix string since it's
        %       just calling a undefined property string
        % 

            testCase.verifyError(@() LtiSystem('DisturbanceMatrixGoneBad', ...
                eye(2)), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectDisturbanceMatrixWrongRows(testCase)
        % Test for failure when calling LtiSystem with mismatch between 
        % disturbance dimension and state dimension
        % 
        % Expected Error: SReachTools:invalidArgs
        % 

            T = 0.5;
            testCase.verifyError(@() LtiSystem('StateMatrix', [1, T; 0, 1], ...
                'DisturbanceMatrix', [T^2], ...
                'DisturbanceSpace', Polyhedron('lb', -1, 'ub', 1)), ...
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectDisturbanceMatrixWrongColumns(testCase)
        % Test for failure when calling LtiSystem with mismatch between 
        % disturbance matrix dimension and disturbance space dimension
        % 
        % Expected Error: SReachTools:invalidArgs
        % 
        % TODO: DisturbanceSpace???
        % 

            T = 0.5;
            testCase.verifyError(@() LtiSystem('StateMatrix', [1, T; 0, 1], ...
                'DisturbanceMatrix', [T^2;T], ...
                'DisturbanceSpace', Polyhedron('lb', [-1;-1], 'ub', [1;1])), ...
                'SReachTools:invalidArgs');
        end
        
        function testIncorrectEmptyDisturbancePolyhedronOneDimDisturbanceMatrix(testCase)
        % Test for failure when calling LtiSystem with disturbance matrix
        % specified but no disturbance
        % 
        % Expected Error: SReachTools:invalidArgs
        % 

            T = 0.5;
            testCase.verifyError(@() LtiSystem('StateMatrix', [1, T; 0, 1], ...
                'DisturbanceMatrix', [T^2;T]), 'SReachTools:invalidArgs');
        end
        
        function testIncorrectRandomVectorBadDim(testCase)
        % Test for failure when calling LtiSystem with mismatch between 
        % disturbance matrix and disturbance
        % 
        % Expected Error: SReachTools:invalidArgs
        % 
        % TODO: Basically the same as the bad column one
        % 

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
        % Test for failure when calling LtiSystem with mistmatch between state
        % dimension and disturbance dimension
        % 
        % Expected Error: SReachTools:invalidArgs
        % 

            T = 0.5;
            testCase.verifyError(@() LtiSystem('StateMatrix', [1, T; 0, 1], ...
                'Disturbance', Polyhedron('lb', -1, 'ub', 1)), ...
                'SReachTools:invalidArgs');
        end
        
        function testCorrectInputPolyhedronOnly(testCase)
        % Test for success of calling LtiSystem with State matrix, input matrix
        % and input space
        % 
        % Expected return: LtiSystem object
        % 

            T = 0.5;
            testCase.verifyClass(LtiSystem('StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2;T], ...
                'InputSpace', Polyhedron('lb', -1, 'ub', 1)), ...
                'LtiSystem');
        end
        
        function testCorrectDisturbancePolyhedronOnly(testCase)
        % Test for success of calling LtiSystem with State matrix, 
        % disurbance matrix and disturbance
        % 
        % Expected return: LtiSystem object
        % 
            T = 0.5;
            testCase.verifyClass(LtiSystem('StateMatrix', [1, T; 0, 1], ...
                'DisturbanceMatrix', [T^2;T], ...
                'Disturbance', Polyhedron('lb', -1, 'ub', 1)), 'LtiSystem');
        end

        function testCorrectArbitrarySystem(testCase)
        % Test for success of calling LtiSystem with State matrix, input matrix,
        % input space, disturbance matrix and disturbance
        % 
        % Expected return: LtiSystem object
        % 

            T = 0.5;
            testCase.verifyClass(LtiSystem('StateMatrix', zeros(2,2), ...
                'InputMatrix', ones(2,4), ...
                'InputSpace', Polyhedron('lb', -ones(4,1), 'ub', ones(4,1)), ...
                'DisturbanceMatrix', ones(2,6), ...
                'Disturbance', Polyhedron('lb', -ones(6,1), 'ub', ones(6,1))), ...
                'LtiSystem');
        end

        function testCorrectDoubleIntegrator(testCase)
        % Test for success of calling LtiSystem with double integrator dynamics
        % 
        % Expected return: LtiSystem object
        % 
        % TODO: Should be changed into the getChainOfInteg testing
        % 

            T = 0.5;
            testCase.verifyClass(LtiSystem('StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2/2;T], ...
                'InputSpace', Polyhedron('lb', -1, 'ub', 1), ...
                'DisturbanceMatrix', [T^2/2;T], ...
                'Disturbance', Polyhedron('lb', -1, 'ub', 1)), 'LtiSystem');
        end
        
        function testCorrectDoubleIntegratorGaussian(testCase)
        % Test for success of calling LtiSystem with double integrator dynamics
        % 
        % Expected return: LtiSystem object
        % 
        % TODO: Should be changed into the getChainOfInteg testing
        % 

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
        % Test for success of calling LtiSystem with double integrator dynamics
        % 
        % Expected return: LtiSystem object
        % 
        % TODO: Should be changed into the getChainOfInteg testing
        % 

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
        
%        %% Test the inherited methods for sanity
%        function testgetHmatMeanCovForXSansInput(testCase)
%            % Input handling tests in LtvSystemTest
%            time_horizon = 10;
%            T = 0.25;
%            umax = 0.75;
%            mean_disturbance = 0;
%            covariance_disturbance = 4;
%            GaussianDisturbance = RandomVector('Gaussian', ...
%                                                         mean_disturbance, ...
%                                                         covariance_disturbance);                                                     
%            sys = LtiSystem(...
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
        % Test method getConcatMats method of LtiSystem
        % 

            % Input handling tests in LtvSystemTest
            time_horizon = 10;
            T = 0.25;
            umax = 0.75;
            dmax = 1;
            mean_disturbance = 0;
            covariance_disturbance = 4;
            GaussianDisturbance = RandomVector('Gaussian', mean_disturbance, ...
                covariance_disturbance);                                                     
            % Load Abar_saved, H_saved, G_matrix_saved expected
            load('./data/getConcatMatsData.mat');
            
            sys = LtiSystem(...
                'StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2; T], ...
                'InputSpace', Polyhedron('lb', -umax, 'ub', umax), ...
                'DisturbanceMatrix', [T^2; T], ...
                'Disturbance', GaussianDisturbance);
            [Z,H,G] = getConcatMats(sys,time_horizon);
            testCase.verifyLessThanOrEqual(sum(sum(abs(Z - Abar_saved))),1e-8);
            testCase.verifyLessThanOrEqual(sum(sum(abs(H - H_matrix_saved))),1e-8);
            testCase.verifyLessThanOrEqual(sum(sum(abs(G - G_matrix_saved))),1e-8);
            % Gave a non-stochastic LtiSystem
            sys = LtiSystem(...
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
            sys = LtiSystem(...
                'StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2; T], ...
                'InputSpace', Polyhedron('lb', -umax, 'ub', umax));
            [Z,H,G] = getConcatMats(sys,time_horizon);
            testCase.verifyLessThanOrEqual(sum(sum(abs(Z - Abar_saved))),1e-8);
            testCase.verifyLessThanOrEqual(sum(sum(abs(H - H_matrix_saved))),1e-8);
            testCase.verifyEqual(size(G),[2*time_horizon 0]);            
            % Control-free
            sys = LtiSystem(...
                'StateMatrix', [1, T; 0, 1], ...
                'DisturbanceMatrix', [T^2; T], ...
                'Disturbance', GaussianDisturbance);
            [Z,H,G] = getConcatMats(sys,time_horizon);
            testCase.verifyLessThanOrEqual(sum(sum(abs(Z - Abar_saved))),1e-8);
            testCase.verifyLessThanOrEqual(sum(sum(abs(G - G_matrix_saved))),1e-8);
            testCase.verifyEqual(size(H),[2*time_horizon 0]);            
        end       
        
        function testgetConcatInputSpace(testCase)
        % Test method getConcatInputSpace method of LtiSystem
        %

            % Input handling tests in LtvSystemTest
            time_horizon = 10;
            umax = 1;
            sys = LtiSystem('StateMatrix', eye(2), ...
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
            sys = getChainOfIntegLtiSystem(2,2, Polyhedron());
            evalc('disp(sys)');
            evalc('disp(sys,''verbose'',true)');
            sys = getChainOfIntegLtiSystem(2,2, Polyhedron(), ...
                RandomVector('gaussian',[0;0],eye(2)));
            evalc('disp(sys)');
            evalc('disp(sys,''verbose'',true)');
            sys = getChainOfIntegLtiSystem(2,2, Polyhedron('lb',-1,'ub',1));
            evalc('disp(sys)');
            evalc('disp(sys,''verbose'',true)');
            sys = getChainOfIntegLtiSystem(2,2, Polyhedron('lb',-1,'ub',1), ...
                RandomVector('gaussian',[0;0],eye(2)));
            evalc('disp(sys)');
            evalc('disp(sys,''verbose'',true)');
        end
    end
end
