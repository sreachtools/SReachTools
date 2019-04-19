classdef SReachDynProgTest < matlab.unittest.TestCase
% Unit tests for SReachDynProg
% ===========================================================================
%
% This function is part of the Stochastic Optimal Control Toolbox.
% License for the use of this function is given in
%      https://sreachtools.github.io/license/
%
%

    methods (Test)
        function testInputHandling(testCase)
            % TODO: Check for a wierd prob_str
            % TODO: Check for non-Gaussian or Ltv System
            % TODO: Check for xinc,uinc
            % TODO: Check for safety_tube (different object, dim)
            % TODO: Check for target_tube (different object, dim, and non-empty)
            % TODO: Check for additional arguments
            % TODO: Check for systems with 4 states
            % TODO: Check for systems with 4 inputs
            % TODO: Check for systems with non-axis aligned hypercuboid
        end

        function dynProgTermBlindTest(testCase)
            % This function blindly runs the dynamic programming --- We will
            % compare the outputs in SReachPointTest.m

            % define the system
            sys = getChainOfIntegLtiSystem(2, 0.1, ...
                    Polyhedron('lb', -0.1, 'ub', 0.1), ...
                    RandomVector('Gaussian', zeros(2,1), 0.01*eye(2)));
            xinc = 0.5;
            uinc = 0.1;
            safe_set = Polyhedron('lb', [-1, -1], 'ub', [1, 1]);
            safety_tube = Tube('viability', safe_set, 5);
            
            % Terminal hitting time problem --- Compared later in
            % chanceConstraint formulation
            SReachDynProg('term', sys, xinc, uinc, safety_tube);

%             % First hitting time problem TODO-first
%             target_set = Polyhedron('lb', 0.5*[-1, -1], 'ub', 0.5*[1, 1]);
%             SReachDynProg('first', sys, xinc, uinc, safety_tube, target_set);
        end
    end
end

