classdef LtiSystemTests < matlab.unittest.TestCase
% Socbox/LtiSystemTests: Unit tests for LtiSystem Class
% ===========================================================================
%
% Unit testing for LtiSystem Class
%
% Usage:
% ------
% tests = LtiSystemTests()
% run(tests)
%   
% ===========================================================================
% 
% LtiSystemTests Methods:
% -----------------------
%   LtiSystemTests/LtiSystemTests - Constructor
%
% ===========================================================================
%
% This function is part of the Stochastic Optimal Control Toolbox.
% License for the use of this function is given in
%      https://github.com/abyvinod/Socbox/blob/master/LICENSE
%
%
    properties (Access = private)
        verbose = false
        T = 0.25
    end

    methods
        function obj = LtiSystemTests(varargin)
        % Socbox/LtiSystemTests/LtiSystemTests: LtiSystemTests constructor
        % ====================================================================
        %
        % Constructor for LtiSystemTests class
        %
        % Usage:
        % ------
        % tests = LtiSystemTests()
        %
        % ====================================================================
        %
        % obj = LtiStystemTests(Name, Value)
        % 
        % Inputs:
        % -------
        %   -----------------------------------------------------
        %   Name              | Value
        %   -----------------------------------------------------
        %   verbose           | true / false
        % 
        % Outputs:
        % --------
        % obj - LtiStystemTests object
        % 
        % ====================================================================
        % 
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://github.com/abyvinod/Socbox/blob/master/LICENSE
        % 
        % 
            p = inputParser();

            p.addParameter('verbose', false, @(x) islogical(x));

            p.parse(varargin{:});

            obj.verbose = p.Results.verbose;
        end
    end

    methods (Test)

        function testInvalidNumOfArgs(test_case)
        % Socbox/LtiSystemTests/testInvalidNumOfArgs: Unit test for invalid 
        % number of arguments
        % =====================================================================
        % 
        % Unit test to test invalid number of arguments by only providing 
        % 'InputMatrix' name with no value. Should throw a 'Socbox:invalidArgs'
        % error.
        % 
        % ====================================================================
        % 
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://github.com/abyvinod/Socbox/blob/master/LICENSE
        % 
        % 

            test_case.verifyError(@() LtiSystem('InputMatrix'), ...
                'Socbox:invalidArgs', ...
                'Did not throw invalidArgs exception');
        end

        function testEmptyLtiSystemError(test_case)
        % Socbox/LtiSystemTests/testEmptyLtiSystemError: Unit test for empty 
        % LtiSystem
        % =====================================================================
        % 
        % Unit test to test empty LtiSystem call; should throw 
        % 'Socbox:invalidArgs' error
        % 
        % ====================================================================
        % 
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://github.com/abyvinod/Socbox/blob/master/LICENSE
        % 
        % 

            test_case.verifyError(@() LtiSystem(), ...
                'Socbox:invalidArgs', ...
                'Did not throw invalidArgs exception');
        end

        function testLtiSystemWithoutStateMatrix(test_case)
        % Socbox/LtiSystemTests/testLtiSystemWithoutStateMatrix: Unit test for
        % LtiSystem call without specified state matrix
        % =====================================================================
        % 
        % Unit test for call to LtiSystem without providing a 'StateMatrix'
        % input; should throw 'Socbox:invalidArgs' error
        % 
        % ====================================================================
        % 
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://github.com/abyvinod/Socbox/blob/master/LICENSE
        % 
        % 

            test_case.verifyError(@() LtiSystem('InputMatrix', eye(2)), ...
                'Socbox:invalidArgs', ...
                'Did not throw invalidArgs exception');
        end

        function testUnhandledArgument(test_case)
        % Socbox/LtiSystemTests/testUnhandledArgument: Unit test for passing
        % and unhandled argument to LtiSystem
        % =====================================================================
        % 
        % Should throw a 'Socbox:invalidArgs' error.
        % 
        % ====================================================================
        % 
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://github.com/abyvinod/Socbox/blob/master/LICENSE
        % 
        % 

            test_case.verifyError(@() LtiSystem('InputMatrixGoneBad', ...
                eye(2)), ...
                'Socbox:invalidArgs', ...
                'Did not throw invalidArgs exception');
        end

        function testInputMatrixRowImbalance(test_case)
        % Socbox/LtiSystemTests/testInputMatrixRowImbalance: Unit test for 
        % imbalanced input matrix
        % =====================================================================
        % 
        % Unit test for call to LtiSystem when the input matrix row count does
        % not match up with state matrix; should throw a 'Socbox:invalidArgs' 
        % error.
        % 
        % =====================================================================
        % 
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://github.com/abyvinod/Socbox/blob/master/LICENSE
        % 
        % 

            T = test_case.T;
            test_case.verifyError(...
                @() LtiSystem('StateMatrix', [1, T; 0, 1], ...
                    'InputMatrix', [T^2], ...
                    'InputSpace', Polyhedron('lb', -1, 'ub', 1)), ...
                'Socbox:invalidArgs', ...
                'Did not throw invalidArgs exception');
        end

        function testInputMatrixColumnImbalance(test_case)
        % Socbox/LtiSystemTests/testInputMatrixColumnImbalance: Unit test for 
        % imbalanced input matrix
        % =====================================================================
        % 
        % Unit test for call to LtiSystem when the input matrix does not have
        % appropriate number of columns, i.e. dim(B*u) ~= dim(x); should throw 
        % a 'Socbox:invalidArgs' error.
        % 
        % =====================================================================
        % 
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://github.com/abyvinod/Socbox/blob/master/LICENSE
        % 
        % 

            T = test_case.T;
            test_case.verifyError(...
                @() LtiSystem('StateMatrix', [1, T; 0, 1], ...
                    'InputMatrix', [T^2;T], ...
                    'InputSpace', Polyhedron('lb', [-1;-1], 'ub', [1;1])), ...
                'Socbox:invalidArgs', ...
                'Did not throw invalidArgs exception'); 
        end

        function testInvalidInputPropertyCombo(test_case)
        % Socbox/LtiSystemTests/testInvalidInputPropertyCombo: Unit test for
        % invalid input proprety combination
        % =====================================================================
        % 
        % Unit test for call to LtiSystem when the input is improperly specified
        % e.g. InputMatrix is provided without and InputSpace; should throw
        % a 'Socbox:invalidArgs' error
        % 
        % =====================================================================
        % 
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://github.com/abyvinod/Socbox/blob/master/LICENSE
        % 
        % 

            T = test_case.T;
            test_case.verifyError(...
                @() LtiSystem('StateMatrix', [1, T; 0, 1], ...
                    'InputMatrix', ones(2,4)), ...
                'Socbox:invalidArgs', ...
                'Did not throw invalidArgs exception');
        end

        function testDisturbanceMatrixRowImbalance(test_case)
        % Socbox/LtiSystemTests/testDisturbanceMatrixRowImbalance: Unit test
        % for disturbance matrix imbalance
        % =====================================================================
        % 
        % Unit test for call to LtiSystem when the disturbance matrix does not
        % have the appropriate number of rows
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

            T = test_case.T;
            test_case.verifyError(...
                @() LtiSystem('StateMatrix', [1, T; 0, 1], ...
                    'DisturbanceMatrix', [T^2], ...
                    'Disturbance', Polyhedron('lb', -1, 'ub', 1)), ...
                'Socbox:invalidArgs', ...
                'Did not throw invalidArgs exception');
        end

        function testDisturbanceMatrixColumnImbalance(test_case)
        % Socbox/LtiSystemTests/testDisturbanceMatrixColumnImbalance: Unit test
        % for disturbance matrix imbalance
        % =====================================================================
        % 
        % Unit test for call to LtiSystem when the disturbance matrix does not
        % have the appropriate number of columns
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

            T = test_case.T;
            test_case.verifyError(...
                @() LtiSystem('StateMatrix', [1, T; 0, 1], ...
                    'DisturbanceMatrix', [T^2;T], ...
                    'Disturbance', Polyhedron('lb', [-1,-1], 'ub', [1,1])), ...
                'Socbox:invalidArgs', ...
                'Did not throw invalidArgs exception');
        end

        function testNonSquareStateMatrix(test_case)
        % Socbox/LtiSystemTests/testNonSquareStateMatrix: Unit test for non
        % square state matrix
        % =====================================================================
        % 
        % Unit test for call to LtiSystem when the state matrix is not square,
        % currently not acceptable because of the need for inversion in the
        % reachability calculations
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

            test_case.verifyError(@() LtiSystem('StateMatrix', [1; 1]), ...
                'Socbox:invalidArgs', ...
                'Did not throw invalidArgs exception');
        end

        function testInvlideDisturbancePropertyCombo(test_case)
        % Socbox/LtiSystemTests/testInvlideDisturbancePropertyCombo: Unit test
        % for improper combination of disturbance properties
        % =====================================================================
        % 
        % Unit test for call to LtiSystem when the combination of disturbance
        % properties provided is invalid, e.g. when a DisturbanceMatrix is 
        % provided without a Disturbance object
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

            T = test_case.T;
            test_case.verifyError(...
                @() LtiSystem('StateMatrix', [1, T; 0, 1], ...
                    'DisturbanceMatrix', ones(2,4)), ...
                'Socbox:invalidArgs', ...
                'Did not throw invalidArgs exception');
        end

        function testStochasticDisturbanceAndMatrixMismatch(test_case)
        % Socbox/LtiSystemTests/testStochasticDisturbanceAndMatrixMismatch: 
        % Unit test for mismatch in Stochastic Disturbance and Matrix
        % =====================================================================
        % 
        % Unit test for call to LtiSystem when the provided dimension of a 
        % StochasticDisturbance object and the Disturbance Matrix to not match
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

            dist_mean = zeros(5, 1);
            dist_cov  = eye(5);
            disturbance = StochasticDisturbance('Gaussian', ...
                                                dist_mean, ...
                                                dist_cov);

            T = test_case.T;
            test_case.verifyError(...
                @() LtiSystem('StateMatrix', [1, T; 0, 1], ...
                    'DisturbanceMatrix', ones(2,4), ...
                    'Disturbance', disturbance), ...
                'Socbox:invalidArgs', ...
                'Did not throw invalidArgs exception');
        end

        function testPartiallyDefinedIntervalDisturbance(test_case)
        % Socbox/LtiSystemTests/testPartiallyDefinedIntervalDisturbance: Unit
        % test for partially (correctly) defined system
        % =====================================================================
        % 
        % Unit test for call to LtiSystem of a partially defined system that is
        % still correct. This test calls a system with a specifed 2x2 
        % StateMatrix and a interval disturbance
        % 
        % =====================================================================
        % 
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://github.com/abyvinod/Socbox/blob/master/LICENSE
        % 
        %

            T = test_case.T;
            test_case.verifyInstanceOf(...
                LtiSystem('StateMatrix', [1, T; 0, 1], ...
                    'Disturbance', Polyhedron('lb', -1, 'ub', [1])), ...
                'LtiSystem');
        end

        function testPartiallyDefinedIntervalInputSpace(test_case)
        % Socbox/LtiSystemTests/testPartiallyDefinedIntervalInputSpace: Unit
        % test for partially (correctly) defined system
        % =====================================================================
        % 
        % Unit test for call to LtiSystem of a partially defined system that is
        % still correct. This test calls a system with a specifed 2x2 
        % StateMatrix and a interval input space
        % 
        % =====================================================================
        % 
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://github.com/abyvinod/Socbox/blob/master/LICENSE
        % 
        %

            T = test_case.T;
            test_case.verifyInstanceOf(...
                LtiSystem('StateMatrix', [1, T; 0, 1], ...
                    'InputSpace', Polyhedron('lb', -1, 'ub', [1])), ...
                'LtiSystem');
        end

        function testCorrectLtiSystemInputOnly(test_case)
        % Socbox/LtiSystemTests/testCorrectLtiSystemInputOnly: Unit test for
        % correct call to LtiSystem
        % =====================================================================
        % 
        % Unit test for complete call to LtiSystem with only state and inputs
        % 
        % =====================================================================
        % 
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://github.com/abyvinod/Socbox/blob/master/LICENSE
        % 
        %
            T = test_case.T;
            test_case.verifyInstanceOf(...
                LtiSystem('StateMatrix', [1, T; 0, 1], ...
                        'InputMatrix', [T^2;T], ...
                        'InputSpace', Polyhedron('lb', -1, 'ub', 1)), ...
                'LtiSystem');
        end

        function testCorrectLtiSystemDisturbanceOnly(test_case)
        % Socbox/LtiSystemTests/testCorrectLtiSystemDisturbanceOnly: Unit test 
        % for correct call to LtiSystem
        % =====================================================================
        % 
        % Unit test for complete call to LtiSystem with only state and 
        % disturbance
        % 
        % =====================================================================
        % 
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://github.com/abyvinod/Socbox/blob/master/LICENSE
        % 
        %

            T = test_case.T;
            test_case.verifyInstanceOf(...
                LtiSystem('StateMatrix', [1, T; 0, 1], ...
                    'DisturbanceMatrix', [T^2;T], ...
                    'Disturbance', Polyhedron('lb', -1, 'ub', 1)), ...
                'LtiSystem');
        end

        function testCorrectLtiSystemArbitrary(test_case)
        % Socbox/LtiSystemTests/testCorrectLtiSystemArbitrary: Unit test for
        % correct call to LtiSystem
        % =====================================================================
        % 
        % Unit test for complete call to LtiSystem with arbitrary (zero and 
        % ones matrices) system definition. All properties provided.
        % 
        % =====================================================================
        % 
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://github.com/abyvinod/Socbox/blob/master/LICENSE
        % 
        %

            sys = LtiSystem('StateMatrix', zeros(2,2), ...
                    'InputMatrix', ones(2,4), ...
                    'InputSpace', Polyhedron('lb', -ones(4,1), ...
                        'ub', ones(4,1)),...
                    'DisturbanceMatrix', ones(2,6), ...
                    'Disturbance', Polyhedron('lb', -ones(6,1), ...
                        'ub', ones(6,1)));

            test_case.verifyInstanceOf(sys, 'LtiSystem');
        end

        function testCorrectLtiSystemDoubleIntegrator(test_case)
        % Socbox/LtiSystemTests/testCorrectLtiSystemDoubleIntegrator: Unit test for
        % correct call to LtiSystem
        % =====================================================================
        % 
        % Unit test for complete call to LtiSystem with double integrator 
        % dynamics
        % 
        % =====================================================================
        % 
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://github.com/abyvinod/Socbox/blob/master/LICENSE
        % 
        %

            T = test_case.T;
            sys = LtiSystem('StateMatrix', [1, T; 0, 1], ...
                    'InputMatrix', [T^2/2;T], ...
                    'InputSpace', Polyhedron('lb', -1, 'ub', 1),...
                    'DisturbanceMatrix', [T^2/2;T], ...
                    'Disturbance', Polyhedron('lb', -1, 'ub', 1));

            test_case.verifyInstanceOf(sys, 'LtiSystem');
        end

        function testCorrestLtiSystemGaussianDisturbance(test_case)
        % Socbox/LtiSystemTests/testCorrestLtiSystemGaussianDisturbance: Unit 
        % test for correct call to LtiSystem with Gaussian Disturbance
        % =====================================================================
        % 
        % Unit test for complete call to LtiSystem with double integrator 
        % dynamics and a Gaussian StochasticDisturbance
        % 
        % =====================================================================
        % 
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://github.com/abyvinod/Socbox/blob/master/LICENSE
        % 
        %

            T = test_case.T;
            dist_mean   = zeros(5,1);
            dist_cov    = eye(5);
            disturbance = StochasticDisturbance('Gaussian',...
                                                dist_mean,...
                                                dist_cov);

            sys = LtiSystem('StateMatrix', [1, T; 0, 1], ...
                            'DisturbanceMatrix', ones(2,5), ...
                            'Disturbance', disturbance);
        end

        function testCorrectDoubleIntWithGaussian(test_case)
        % Socbox/LtiSystemTests/testCorrestLtiSystemGaussianDisturbance: Unit 
        % test for correct call to LtiSystem with Gaussian Disturbance
        % =====================================================================
        % 
        % Unit test for complete call to LtiSystem with double integrator 
        % dynamics and a Gaussian StochasticDisturbance
        % 
        % =====================================================================
        % 
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://github.com/abyvinod/Socbox/blob/master/LICENSE
        % 
        %
        
            T = test_case.T;
            dist_mean   = 0;
            dist_cov    = 1;
            disturbance = StochasticDisturbance('Gaussian',...
                                                dist_mean,...
                                                dist_cov);

            sys = LtiSystem('StateMatrix', [1, T; 0, 1], ...
                            'InputMatrix', [T^2/2;T], ...
                            'InputSpace', Polyhedron('lb', -1, 'ub', 1),...
                            'DisturbanceMatrix', [T^2/2;T], ...
                            'Disturbance', disturbance);
        end
    end
end