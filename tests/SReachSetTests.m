classdef SReachSetTests < matlab.unittest.TestCase
    methods (Test)
         
        function testPointBasedSetComputations(test_case)
            [sys, safety_tube] = test_case.getDI();
            theta = linspace(0,2*pi,3);
            options = SReachSetOptions('term','chance-open', ...
                'set_of_dir_vecs',[cos(theta);sin(theta)], ...
                'init_safe_set_affine',Polyhedron(), 'verbose',0);
            optionsGp = SReachSetOptions('term','genzps-open', ...
                'set_of_dir_vecs',[cos(theta);sin(theta)], ...
                'init_safe_set_affine',Polyhedron(), 'verbose', 0, ...
                'desired_accuracy', 0.1);

            %% Non-empty case - chance-open
            [level_set, extra_info] = SReachSet('term','chance-open', sys, ...
                0.8, safety_tube, options);
            test_case.verifyInstanceOf(level_set, 'Polyhedron');
            % extra_info
            
            %% Non-empty case - genzps-open
            [level_set, extra_info] = SReachSet('term','genzps-open', sys, ...
                0.8, safety_tube, optionsGp);
            test_case.verifyInstanceOf(level_set, 'Polyhedron');             
        end
    
        function testLagrangianUnderapproxAndSReachSetLagBset(test_case)
            [sys, safety_tube] = test_case.getDI();

            % Test underapproximation using polytope option
            bound_vec = ones(sys.dist_dim,1);
            luOpts_poly = SReachSetOptions('term', 'lag-under', ...
                'bound_set_method', 'polytope', 'template_polytope',...
                Polyhedron('lb',-bound_vec,'ub',bound_vec));
            luSet = SReachSet('term', 'lag-under', sys, 0.8, safety_tube, ...
                luOpts_poly);
            test_case.verifyInstanceOf(luSet, 'Polyhedron');

            % Test underapproximation using ellipsoid option
            n_dim = (sys.state_dim + sys.input_dim);
            luOpts_ellipsoid = SReachSetOptions('term', 'lag-under', ...
                'bound_set_method', 'ellipsoid', ...
                'system', sys, 'verbose', 0, ...
                'n_vertices',2^n_dim * 10 + 2*n_dim);
            luSet = SReachSet('term', 'lag-under', sys, 0.8, safety_tube, ...
                luOpts_ellipsoid);

            test_case.verifyInstanceOf(luSet, 'Polyhedron');
        end 
        
        function testLagrangianOverapprox(test_case)
            [sys, safety_tube] = test_case.getDI();

            % Test overapproximation using polytope option
            bound_vec = ones(sys.dist_dim,1);
            loOpts = SReachSetOptions('term', 'lag-over', ...
                'bound_set_method', 'polytope', 'template_polytope',...
                Polyhedron('lb',-bound_vec,'ub',bound_vec));
            loSet = SReachSet('term', 'lag-over', sys, 0.8, safety_tube, ...
                loOpts);

            test_case.verifyInstanceOf(loSet, 'Polyhedron');

        end 
        
        function getBsetWithProbTest(test_case)
            [sys, safety_tube] = test_case.getDI();
            dist = 100*sys.dist;
            
            bound_vec = ones(sys.dist_dim,1);
            temp_polytope = Polyhedron('lb',-bound_vec,'ub',bound_vec);
            onestep_prob_thresh = 0.8^(1/(length(safety_tube)-1));
            
            % Test getBsetWithProb explicitly
            desired_accuracy = 1e-2;
            bounded_set = getBsetWithProb(dist, temp_polytope,...
                onestep_prob_thresh, desired_accuracy, 0);
            n_particles_test = 1e6;
            mcarlo_sims = dist.getRealizations(n_particles_test);
            count_contains = bounded_set.contains(mcarlo_sims);
            prob_test = sum(count_contains)/n_particles_test;
            % The polytope returned must have a probability no smaller than
            % the requested threshold
            test_case.verifyLessThanOrEqual(onestep_prob_thresh,prob_test);          
        end                        
        function testSReachSetGpOEmpty(test_case)
            T = 0.25;

            % Safe set K |x1| < 1, |x2| < 1
            safe_set = Polyhedron('lb', [-1; -1], 'ub', [1; 1]);
            time_horizon = 6;
            % target_tube = {K, K, K, K, K, K};
            safety_tube = Tube('viability', safe_set, time_horizon);

            % Input Space
            U = Polyhedron('lb', -0.1, 'ub', 0.1);

            %% Empty case: 0.999 => emptyset, but non-trivial Wmax
            % getDI only had 1e-3 covariance
            sysEmpty = LtiSystem('StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2/2; T], ...
                'InputSpace', U, ...
                'DisturbanceMatrix', eye(2), ...
                'Disturbance', RandomVector('Gaussian', zeros(2,1), ...
                    10*eye(2)));
                
            % Define options
            theta = linspace(0,2*pi,3);            
            optionsGp = SReachSetOptions('term','genzps-open', ...
                'set_of_dir_vecs',[cos(theta);sin(theta)], ...
                'init_safe_set_affine',Polyhedron(), 'verbose', 0, ...
                'desired_accuracy', 0.1);

            %% Empty case - genzps-open           
            [level_set, extra_info] = SReachSet('term','genzps-open', ...
                sysEmpty, 0.999, safety_tube, optionsGp);
            test_case.verifyInstanceOf(level_set, 'Polyhedron');
            test_case.verifyTrue(level_set.isEmptySet(), ...
                'Expected empty stochastic reach set');
            test_case.verifyEqual(length(extra_info.extra_info_cco), 2);
            test_case.verifyTrue(abs(extra_info(1).xmax_reach_prob-0.1) <0.1,...
                'Expected extra_info(1).xmax_reach_prob to be close to 0.1');
        end
        function testSReachSetCcOAllComputeStylesAndEmpty(test_case)

            T = 0.25;

            % Safe set K |x1| < 1, |x2| < 1
            safe_set = Polyhedron('lb', [-1; -1], 'ub', [1; 1]);
            time_horizon = 6;
            % target_tube = {K, K, K, K, K, K};
            safety_tube = Tube('viability', safe_set, time_horizon);

            % Input Space
            U = Polyhedron('lb', -0.1, 'ub', 0.1);

            %% Empty case: 0.999 => emptyset, but non-trivial Wmax
            % getDI only had 1e-3 covariance
            sys = LtiSystem('StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2/2; T], ...
                'InputSpace', U, ...
                'DisturbanceMatrix', eye(2), ...
                'Disturbance', RandomVector('Gaussian', zeros(2,1), ...
                    1e-2*eye(2)));

            theta = linspace(0,2*pi,3);            
            
            %% Empty case - chance-open
            options = SReachSetOptions('term','chance-open', ...
                'set_of_dir_vecs',[cos(theta);sin(theta)], ...
                'init_safe_set_affine',Polyhedron(), 'verbose', 0);
            [level_set, extra_info] = SReachSet('term','chance-open', sys, ...
                0.999, safety_tube, options);
            test_case.verifyInstanceOf(level_set, 'Polyhedron');
            test_case.verifyTrue(level_set.isEmptySet(), ...
                'Expected empty stochastic reach set');
            test_case.verifyEqual(length(extra_info), 2);
            test_case.verifyTrue(extra_info(1).xmax_reach_prob < 0.999, ...
                'Expected non-trivial extra_info(1).xmax_reach_prob');

            %% Non-empty case - chance-open (wmax only)
            options = SReachSetOptions('term','chance-open', ...
                'set_of_dir_vecs',[cos(theta);sin(theta)], ...
                'init_safe_set_affine',Polyhedron(), 'verbose', 0, ...
                'compute_style','max_safe_init');
            [level_set, extra_info] = SReachSet('term','chance-open', sys, ...
                0.8, safety_tube, options);
            test_case.verifyInstanceOf(level_set, 'Polyhedron');
            test_case.verifyTrue(~level_set.isEmptySet(), ...
                'Expected non-empty stochastic reach set');
            test_case.verifyEqual(length(extra_info), 2);
            test_case.verifyTrue(extra_info(1).xmax_reach_prob > 0.8, ...
                'Expected non-trivial extra_info(1).xmax_reach_prob');

            %% Non-empty case - chance-open (cheby only)
            options = SReachSetOptions('term','chance-open', ...
                'set_of_dir_vecs',[cos(theta);sin(theta)], ...
                'init_safe_set_affine',Polyhedron(), 'verbose', 0, ...
                'compute_style','cheby');
            [level_set, extra_info] = SReachSet('term','chance-open', sys, ...
                0.8, safety_tube, options);
            test_case.verifyInstanceOf(level_set, 'Polyhedron');
            test_case.verifyTrue(~level_set.isEmptySet(), ...
                'Expected non-empty stochastic reach set');
            test_case.verifyEqual(length(extra_info), 2);
            test_case.verifyTrue(extra_info(1).xmax_reach_prob > 0.8, ...
                'Expected non-trivial extra_info(1).xmax_reach_prob');

            %% Non-empty case - chance-open (all)
            options = SReachSetOptions('term','chance-open', ...
                'set_of_dir_vecs',[cos(theta);sin(theta)], ...
                'init_safe_set_affine',Polyhedron(), 'verbose', 0, ...
                'compute_style','all');
            [level_set, extra_info] = SReachSet('term','chance-open', sys, ...
                0.8, safety_tube, options);
            test_case.verifyInstanceOf(level_set, 'Polyhedron');
            test_case.verifyTrue(~level_set.isEmptySet(), ...
                'Expected non-empty stochastic reach set');
            test_case.verifyEqual(length(extra_info), 2);
            
            %% Empty case: 0.999 => emptyset for high stochastic system and Wmax
            %% is trival
            sysEmpty = LtiSystem('StateMatrix', sys.state_mat, ...
                                 'InputMatrix', sys.input_mat, ...
                                 'InputSpace', sys.input_space, ...
                                 'DisturbanceMatrix', sys.dist_mat, ...
                                 'Disturbance', RandomVector('Gaussian', ...
                                    zeros(2,1), 10*eye(2)));

            %% Empty case - chance-open
            options = SReachSetOptions('term','chance-open', ...
                'set_of_dir_vecs',[cos(theta);sin(theta)], ...
                'init_safe_set_affine',Polyhedron(), 'verbose', 0);
            [level_set, extra_info] = SReachSet(...
                'term','chance-open', sysEmpty, 0.999, safety_tube, options);
            test_case.verifyInstanceOf(level_set, 'Polyhedron');
            test_case.verifyTrue(level_set.isEmptySet(), ...
                'Expected empty stochastic reach set');
            test_case.verifyEqual(length(extra_info), 2);
            test_case.verifyTrue(isempty(extra_info(1).xmax_reach_prob), ...
                'Expected extra_info(1).xmax_reach_prob to be empty');
        end
    end
    methods (Static)
        function [sys, safety_tube] = getDI()
            T = 0.25;

            % Safe set K |x1| < 1, |x2| < 1
            safe_set = Polyhedron('lb', [-1; -1], 'ub', [1; 1]);
            time_horizon = 6;

            % Input Space
            U = Polyhedron('lb', -0.1, 'ub', 0.1);

            sys = LtiSystem('StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2/2; T], ...
                'InputSpace', U, ...
                'DisturbanceMatrix', eye(2), ...
                'Disturbance', RandomVector('Gaussian', zeros(2,1), ...
                    1e-3*eye(2)));

            % target_tube = {K, K, K, K, K, K};
            safety_tube = Tube('viability', safe_set, time_horizon);
        end
        
        function [sysCwh, safety_tube] = getCwh()
            sysCwh = getCwhLtiSystem(4, ...
                        Polyhedron('lb',-0.01*[1;1],'ub', 0.01*[1;1]), ...
                        RandomVector('Gaussian', zeros(4,1), ...
                            diag([1e-4, 1e-4, 5e-8, 5e-8])));
            %% Safe set definition --- LoS cone |x|<=y and y\in[0,ymax] and 
            % |vx|<=vxmax and 
            %% |vy|<=vymax
            time_horizon = 5;
            ymax = 2;
            vxmax = 0.5;
            vymax = 0.5;
            A_safe_set = [1, 1, 0, 0;           
                         -1, 1, 0, 0; 
                          0, -1, 0, 0;
                          0, 0, 1,0;
                          0, 0,-1,0;
                          0, 0, 0,1;
                          0, 0, 0,-1];
            b_safe_set = [0;
                          0;
                          ymax;
                          vxmax;
                          vxmax;
                          vymax;
                          vymax];
            safe_set = Polyhedron(A_safe_set, b_safe_set);

            %% Target set --- Box [-0.1,0.1]x[-0.1,0]x[-0.01,0.01]x[-0.01,0.01]
            target_set = Polyhedron('lb', [-0.1; -0.1; -0.01; -0.01], ...
                                    'ub', [0.1; 0; 0.01; 0.01]);
            safety_tube = Tube('reach-avoid',safe_set, target_set, ...
                time_horizon);                    
        end
    end                
end
