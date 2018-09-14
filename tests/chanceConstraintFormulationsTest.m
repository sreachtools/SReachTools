classdef chanceConstraintFormulationsTest < matlab.unittest.TestCase

    methods (Test)
        function testChanceConstraintFormulation(testCase)
            
            n_test_pts = 10;
            
            dyn_prog_xinc = 0.05;
            dyn_prog_uinc = 0.5;
            
            time_horizon = 5;
            
            n_intg = 2;
            umax = 1;
            xmax = [1,1];
            sampling_time = 0.1;
            dist_cov = 0.001;
            sys = getChainOfIntegLtiSystem(n_intg,...
                sampling_time,...
                Polyhedron('lb',-umax,'ub',umax),...
                RandomVector('Gaussian', zeros(n_intg,1), dist_cov * eye(n_intg)));

            %% Setup the target tube
            % safe set definition
            safe_set = Polyhedron('lb', -xmax, 'ub', xmax);
            % target tube definition
            target_tube = TargetTube('viability', safe_set, time_horizon);

            [prob_x, ~, grid_x] = getDynProgSolForTargetTube(sys, dyn_prog_xinc, dyn_prog_uinc, target_tube);
            
            % No. of grid points
            n_grid_x = length(grid_x);

            % Randomly generate n_test_pts from the middle of the grid
            rand_indices = datasample(find(max(abs(grid_x'))<xmax(1)/3 == 1),n_test_pts,'Replace',false);
            true_probability = prob_x(rand_indices);
            
%             lb_safe_prob = zeros(n_test_pts,1);
%             for indx=1:n_test_pts
%                 initial_state = grid_x(indx,:)';
%                 lb_safe_prob(indx) = getLowerBoundStochReachAvoid(sys,...
%                                              initial_state,...
%                                              target_tube,...
%                                              'cccpwl');
%             end
%             testCase.verifyTrue(all(lb_safe_prob<=true_probability),'Not a lower bound?');
%             testCase.verifyTrue(all(lb_safe_prob>=0),'Shouldn''t be infeasible!');
        end
    end
end
