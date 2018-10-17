classdef BoundedDisturbanceTests < matlab.unittest.TestCase
% SReachTools/BoundedDisturbanceTests: Unit tests for bounded disturbances
% ===========================================================================
%
% Unit testing for bounded disturbances
%
% Usage:
% ------
% tests = BoundedDisturbanceTests()
% run(tests)
%
% ===========================================================================
%
% This function is part of the Stochastic Optimal Control Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
%
%

    methods (Test)
        function testBoundedEllipseRandomDirections(test_case)
        % SReachTools/BoundedDisturbanceTests/testBoundedEllipseRandomDirections: 
        % Unit test for getBoundedSetForDisturbance for random ellipsoid 
        % generation
        % ====================================================================
        %
        % Unit test for getBoundedSetForDisturbance using ellipsoid 
        % overapproximation via random direction selection
        %   
        % ====================================================================
        %
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        %
        %   

            T = 0.25;
            umax = 0.75;
            dmax = 0.097;

            sys = LtiSystem(...
                'StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2; T], ...
                'InputSpace', Polyhedron('lb', -umax, 'ub', umax), ...
                'DisturbanceMatrix', eye(2), ...
                'Disturbance', RandomVector('Gaussian', zeros(2,1),...
                    0.01* eye(2)));

            beta = 0.7;
            time_horizon = 5;

            bounded_dist = getBoundedSetForDisturbance(...
                sys.dist, ...
                time_horizon, ...
                beta, ...
                'random', ...
                100);

            test_case.verifyInstanceOf(bounded_dist, 'Polyhedron');
        end

        function testLoadBoundedDisturbance(test_case)
        % SReachTools/BoundedDisturbanceTests/testBoundedEllipseRandomDirections: 
        % Unit test for getBoundedSetForDisturbance for random ellipsoid 
        % generation
        % ====================================================================
        %
        % Unit test for getBoundedSetForDisturbance using load from file method
        %   
        % ====================================================================
        %
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        %
        %

            bounded_dist = getBoundedSetForDisturbance(...
                [],...
                [],...
                [],... 
                'load', ...
                'data/cwhUnderapproxBoundeSet.mat');

            test_case.verifyInstanceOf(bounded_dist, 'Polyhedron');
        end
    end
end
