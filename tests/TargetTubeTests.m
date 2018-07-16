classdef TargetTubeTests < matlab.unittest.TestCase
% Socbox/TargetTubeTests: Unit tests for bounded disturbances
% ===========================================================================
%
% Unit testing for bounded disturbances
%
% Usage:
% ------
% tests = TargetTubeTests()
% run(tests)
%
% ===========================================================================
%
% This function is part of the Stochastic Optimal Control Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/Socbox/blob/master/LICENSE
%
%

    methods (Test)
        function testViabilityTargetTube(test_case)
            time_horizon = 5;
            tt = TargetTube('viability', Polyhedron('lb', [-1, -1], ...
                'ub', [1, 1]), time_horizon);

            test_case.verifyInstanceOf(tt, 'TargetTube');
            test_case.verifyLength(tt, time_horizon+1);
        end

        function testReachAvoidTargetTube(test_case)
            time_horizon = 5;

            target_set = Polyhedron('lb', [-0.5, -0.5], 'ub', [0.5, 0.5]);
            safe_set = Polyhedron('lb', [-1, -1], 'ub', [1, 1]);

            tt = TargetTube('reach-avoid', safe_set, target_set, time_horizon);

            test_case.verifyInstanceOf(tt, 'TargetTube');
            test_case.verifyLength(tt, time_horizon+1);
        end

        function testInvalidTypeCall(test_case)
            time_horizon = 5;

            target_set = Polyhedron('lb', [-0.5, -0.5], 'ub', [0.5, 0.5]);
            safe_set = Polyhedron('lb', [-1, -1], 'ub', [1, 1]);

            test_case.verifyError( ...
                @() TargetTube('invalid', safe_set, target_set, ...
                    time_horizon), ...
                ?SrtInvalidArgsError);
        end

        function testGenericTargetTube(test_case)
            time_horizon = 5;

            target_sets = cell(1, time_horizon);
            for lv = 1:time_horizon
                target_sets{lv} = Polyhedron('lb', -rand(), 'ub', rand());
            end

            tt = TargetTube(target_sets{:});

            test_case.verifyInstanceOf(tt, 'TargetTube');
            test_case.verifyLength(tt, time_horizon);
        end

        function testNonpolyFirstInput(test_case)
            test_case.verifyError(@() TargetTube(3), 'MATLAB:invalidType');
        end

        function testDifferentDimensionPolyhedron(test_case)
            test_case.verifyError( ...
                @() TargetTube( ...
                    Polyhedron(eye(2), ones(2, 1)), ...
                    Polyhedron(1, 1)), ...
                ?SrtInvalidArgsError);
        end
    end
end
