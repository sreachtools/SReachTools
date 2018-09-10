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
        
        function testConcatTube(testCase)
            time_horizon = 10;
            xmax = 1;
            
            % Target tubes has polyhedra T_0, T_1, ..., T_{time_horizon}
            target_tube = TargetTube('viability',Polyhedron('lb',-xmax,'ub',xmax), time_horizon);
            
            % concat using all the time steps --- empty
            [concat_target_tube_A, concat_target_tube_b] = target_tube.concat();            
            obtained_polyhedron = Polyhedron('H',[concat_target_tube_A, ...
                            concat_target_tube_b]);
            expected_polyhedron = Polyhedron('lb', -xmax * ones(time_horizon+1,1), ...
                       'ub',  xmax * ones(time_horizon+1,1));
            testCase.verifyTrue(obtained_polyhedron == expected_polyhedron);
            
            % concat using all the time steps --- both ends
            [concat_target_tube_A, concat_target_tube_b] = target_tube.concat([1 time_horizon+1]);            
            obtained_polyhedron = Polyhedron('H',[concat_target_tube_A, ...
                            concat_target_tube_b]);
            expected_polyhedron = Polyhedron('lb', -xmax * ones(time_horizon+1,1), ...
                       'ub',  xmax * ones(time_horizon+1,1));
            testCase.verifyTrue(obtained_polyhedron == expected_polyhedron);
            
            % concat using 1 to time horizon
            xmax = [1 10 2 3];
            target_tube = TargetTube(Polyhedron('lb',-xmax(1),'ub',xmax(1)), ...
                                     Polyhedron('lb',-xmax(2),'ub',xmax(2)), ...
                                     Polyhedron('lb',-xmax(3),'ub',xmax(3)), ...
                                     Polyhedron('lb',-xmax(4),'ub',xmax(4)));
            [concat_target_tube_A, concat_target_tube_b] = target_tube.concat([4 4]);            
            obtained_polyhedron = Polyhedron('H',[concat_target_tube_A, ...
                            concat_target_tube_b]);
            expected_polyhedron = Polyhedron('lb', -xmax(4), 'ub',  xmax(4));
            testCase.verifyTrue(obtained_polyhedron == expected_polyhedron);
            
            % concat using 1 to time horizon
            xmax = [1 10 2 3];
            target_tube = TargetTube(Polyhedron('lb',-xmax(1),'ub',xmax(1)), ...
                                     Polyhedron('lb',-xmax(2),'ub',xmax(2)), ...
                                     Polyhedron('lb',-xmax(3),'ub',xmax(3)), ...
                                     Polyhedron('lb',-xmax(4),'ub',xmax(4)));
            [concat_target_tube_A, concat_target_tube_b] = target_tube.concat([2 4]);            
            obtained_polyhedron = Polyhedron('H',[concat_target_tube_A, ...
                            concat_target_tube_b]);
            expected_polyhedron = Polyhedron('lb', -xmax(2:4), 'ub',  xmax(2:4));
            testCase.verifyTrue(obtained_polyhedron == expected_polyhedron);
            
            % concat using 1 to time_horizon + 1 (will throw error)
            testCase.verifyError(@() target_tube.concat([1 5]),'SReachTools:invalidArgs');            
            testCase.verifyError(@() target_tube.concat([-1 5]),'SReachTools:invalidArgs');            
            testCase.verifyError(@() target_tube.concat([-1 2 5]),'SReachTools:invalidArgs');            
        end
        
        function testContains(testCase)
            time_horizon = 10;
            xmax = 1;
            
            % Target tubes has polyhedra T_0, T_1, ..., T_{time_horizon}
            target_tube = TargetTube('viability',Polyhedron('lb',-xmax,'ub',xmax), time_horizon);
            
            %% Checking input handling
            % Empty vector
            testCase.verifyError(@() target_tube.contains(repmat(0,time_horizon+1,0)),'SReachTools:invalidArgs');            
            % Incorrect time horizon
            testCase.verifyError(@() target_tube.contains(repmat(0,1,2)),'SReachTools:invalidArgs');            
            %% Checking correct behaviour
            % Single concatenated state trajectories
            testCase.verifyTrue(target_tube.contains(repmat(0,time_horizon+1,1)));
            % Multiple concatenated state trajectories
            testCase.verifyEqual(target_tube.contains(repmat([xmax/2 xmax],time_horizon+1,1)),logical([1 1]));
            % Single concatenated state trajectory outside
            testCase.verifyFalse(target_tube.contains(repmat(xmax*2,time_horizon+1,1)));
            % Multiple concatenated state trajectories with one of them outside
            testCase.verifyEqual(target_tube.contains(repmat([xmax/2 xmax 2*xmax],time_horizon+1,1)),logical([1 1 0]));            
        end
    end
end
