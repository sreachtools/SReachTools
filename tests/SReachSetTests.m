classdef SReachSetTests < matlab.unittest.TestCase
    methods (Test)
%         function testLagrangian(test_case)
%             [sys, safety_tube] = test_case.getDI();
%             %% Random option
%             % Underaproximation
%             options = SReachSetOptions('term','lag-under',...
%                 'bound_set_method', 'random',...
%                 'num_dirs',100);
%             level_set = SReachSet('term','lag-under', sys, 0.8, safety_tube,...
%                 options);
%             test_case.verifyInstanceOf(level_set, 'Polyhedron');
%             
%             % Overaproximation
%             options = SReachSetOptions('term','lag-over',...
%                 'bound_set_method', 'random',...
%                 'num_dirs',100);
%             level_set = SReachSet('term','lag-over', sys, 0.8, safety_tube,...
%                 options);
%             test_case.verifyInstanceOf(level_set, 'Polyhedron');
%             
%             %% TODO: Repeat the same for box, optim-box options
%             [sys_cwh, safety_tube_cwh] = test_case.getCwh();
%             %% Load option
%             % Underaproximation
%             options = SReachSetOptions('term','lag-under',...
%                 'bound_set_method', 'load',...
%                 'load_str','data/cwhUnderapproxBoundeSet.mat');
%             level_set = SReachSet('term','lag-under', sys_cwh, 0.8,...
%                 safety_tube_cwh, options);
%             test_case.verifyInstanceOf(level_set, 'Polyhedron');
%             
% %             % Overaproximation
%             options = SReachSetOptions('term','lag-over',...
%                 'bound_set_method', 'load',...
%                 'load_str','data/cwhUnderapproxBoundeSet.mat');
%             level_set = SReachSet('term','lag-over', sys_cwh, 0.8,...
%                 safety_tube_cwh, options);
%             test_case.verifyInstanceOf(level_set, 'Polyhedron');  
% 
%             [sys, safety_tube] = test_case.getDI();
%             %% Random option
%             % Underaproximation
%             options = SReachSetOptions('term','lag-under',...
%                 'bound_set_method', 'box',...
%                 'num_dirs',100);
%             level_set = SReachSet('term','lag-under', sys, 0.8, safety_tube,...
%                 options);
%             test_case.verifyInstanceOf(level_set, 'Polyhedron');
%             
%             % Overaproximation
%             options = SReachSetOptions('term','lag-over',...
%                 'bound_set_method', 'box',...
%                 'num_dirs',100);
%             level_set = SReachSet('term','lag-over', sys, 0.8, safety_tube,...
%                 options);
%             test_case.verifyInstanceOf(level_set, 'Polyhedron');                      
%         end
         
        function testPointBased(test_case)
            [sys, safety_tube] = test_case.getDI();
            theta = linspace(0,2*pi,11);
            options = SReachSetOptions('term','chance-open',...
                'set_of_dir_vecs',[cos(theta);sin(theta)],...
                'init_safe_set_affine',Polyhedron(), 'verbose',0);
            % Non-empty case
            level_set = SReachSet('term','chance-open', sys, 0.8,...
                safety_tube, options);
            test_case.verifyInstanceOf(level_set, 'Polyhedron');
            [level_set, extra_info] = SReachSet(...
                'term','chance-open', sys, 0.8, safety_tube, options);

            % Empty case: 0.999 => emptyset for high stochastic system
            sysEmpty = LtiSystem('StateMatrix', sys.state_mat, ...
                                 'InputMatrix', sys.input_mat, ...
                                 'InputSpace', sys.input_space, ...
                                 'DisturbanceMatrix', sys.dist_mat, ...
                                 'Disturbance', RandomVector('Gaussian',...
                                    zeros(2,1), 10*eye(2)));
            level_set = SReachSet('term','chance-open', sysEmpty, 0.999,...
                safety_tube, options);
            test_case.verifyEqual(level_set, Polyhedron.emptySet(2),...
                'Should be empty');
            [level_set, extra_info] = SReachSet(...
                'term','chance-open', sysEmpty, 0.999, safety_tube, options);
            test_case.verifyEqual(level_set, Polyhedron.emptySet(2),...
                'Should be empty');
            test_case.verifyLessThanOrEqual(...
                extra_info(1).xmax_reach_prob,0.8,'Low max safety');
            test_case.verifyEqual(length(extra_info),1,'Empty cheby');
%             figure();
%             plot(safety_tube(1),'alpha',0.3);
%             hold on
%             plot(level_set);
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
                'Disturbance', RandomVector('Gaussian', zeros(2,1), 5e-3*eye(2)));

            % target_tube = {K, K, K, K, K, K};
            safety_tube = Tube('viability', safe_set, time_horizon);
        end
        
        function [sysCwh, safety_tube] = getCwh()
            sysCwh = getCwhLtiSystem(4,...
                        Polyhedron('lb',-0.01*[1;1],'ub', 0.01*[1;1]), ...
                        RandomVector('Gaussian', zeros(4,1), ...
                            diag([1e-4, 1e-4, 5e-8, 5e-8])));
            %% Safe set definition --- LoS cone |x|<=y and y\in[0,ymax] and |vx|<=vxmax and 
            %% |vy|<=vymax
            time_horizon = 5;
            ymax=2;
            vxmax=0.5;
            vymax=0.5;
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
            target_set = Polyhedron('lb', [-0.1; -0.1; -0.01; -0.01],...
                                    'ub', [0.1; 0; 0.01; 0.01]);
            safety_tube = Tube('reach-avoid',safe_set, target_set,...
                time_horizon);                    
        end
    end            
end
