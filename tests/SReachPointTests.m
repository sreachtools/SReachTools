classdef SReachPointTests < matlab.unittest.TestCase
    methods (Test)
         
        function testPointChanceOpen(test_case)
            [sys, safety_tube, initial_state] = test_case.getDI();
            
            [lb_stoch_prob, input_vec, risk_state, extra_info] = SReachPoint(...
                'term','chance-open', sys, initial_state, safety_tube);
        end 
        
        function testPointParticleOpen(test_case)
            [sys, safety_tube, initial_state] = test_case.getDI();
            
            test_case.verifyError(@() SReachPointOptions(...
                    'term','particle-open','n_particles',300),...
                    'SReachTools:invalidArgs');
            % Use default options
            [lb_stoch_prob, input_vec] = SReachPoint('term','particle-open', ...
                sys, initial_state, safety_tube);
        end 
        
        function testPointGenzpsOpen(test_case)
            [sys, safety_tube, initial_state] = test_case.getDI();
            
            [lb_stoch_prob, input_vec] = SReachPoint('term','genzps-open', ...
                sys, initial_state, safety_tube);
        end 
        
        function testPointVoronoiOpen(test_case)
            [sys, safety_tube, initial_state] = test_case.getDI();
            
            options_bad = SReachPointOptions('term','voronoi-open',...
                    'max_overapprox_err', 1e-3);
            test_case.verifyError(@() SReachPoint('term','voronoi-open', sys,...
                initial_state, safety_tube, options_bad),...
                'SReachTools:invalidArgs');
            
            % Working three option
            [lb_stoch_prob, input_vec, extra_info] = SReachPoint('term', ...
                'voronoi-open', sys,initial_state, safety_tube);            
        end 
        
        function testPointChanceAffine(test_case)
            [sys, safety_tube, initial_state] = test_case.getDI();
            
            % Forgotten max_input_viol_prob for chance-affine
            test_case.verifyError(@() SReachPointOptions(...
                'term','chance-affine'),'SReachTools:invalidArgs');
            test_case.verifyError(@() SReachPoint('term','chance-affine', ...
                sys, initial_state, safety_tube),'SReachTools:invalidArgs');
            
            % Working example
            options = SReachPointOptions('term','chance-affine', ...
                'max_input_viol_prob', 2e-1, 'verbose', 0);
            [lb_stoch_prob, input_vec, input_gain, risk_state, risk_input] = ...
                SReachPoint('term','chance-affine', sys, initial_state, ...
                safety_tube, options);
            
            sys_bad = LtiSystem('StateMatrix', sys.state_mat, 'InputMatrix', ...
                sys.input_mat, 'InputSpace', sys.input_space, 'Disturbance', ...
                sys.dist, 'DisturbanceMatrix', magic(sys.dist_dim));
            test_case.verifyError(@() SReachPoint('term','chance-affine', ...
                sys_bad, initial_state, safety_tube),'SReachTools:invalidArgs');

        end  

        function testPointChanceAffineUni(test_case)
            [sys, safety_tube, initial_state] = test_case.getDI();
            
            % Forgotten max_input_viol_prob for chance-affine
            test_case.verifyError(@() SReachPointOptions(...
                'term','chance-affine-uni'),'SReachTools:invalidArgs');
            test_case.verifyError(@() SReachPoint('term','chance-affine-uni',...
                sys, initial_state, safety_tube),'SReachTools:invalidArgs');
            
            % Working example --- Had to set to [0;0] to avoid a trivial bound
            options = SReachPointOptions('term','chance-affine-uni', ...
                'max_input_viol_prob', 2e-1, 'verbose', 0);
            [lb_stoch_prob, input_vec, input_gain, risk_state, risk_input] = ...
                SReachPoint('term','chance-affine-uni', sys, [0;0], ...
                safety_tube, options);
            
            sys_bad = LtiSystem('StateMatrix', sys.state_mat, 'InputMatrix', ...
                sys.input_mat, 'InputSpace', sys.input_space, 'Disturbance', ...
                sys.dist, 'DisturbanceMatrix', magic(sys.dist_dim));
            test_case.verifyError(@() SReachPoint('term','chance-affine-uni', ...
                sys_bad, initial_state, safety_tube),'SReachTools:invalidArgs');

        end  

%         function testPointVoronoiAffine(test_case)
%             [sys, safety_tube, initial_state] = test_case.getDI();
%             
%             % Forgotten max_input_viol_prob for chance-affine
%             test_case.verifyError(@() SReachPointOptions(...
%                 'term','voronoi-affine'),'SReachTools:invalidArgs');
%             test_case.verifyError(@() SReachPoint('term','voronoi-affine', ...
%                 sys, initial_state, safety_tube),'SReachTools:invalidArgs');
%             
%             % Bad constraint: TODO 1 - Du + d<= 1 - d
%             
%             % Working example
%             options = SReachPointOptions('term','voronoi-affine',...
%                 'verbose', 0, 'max_input_viol_prob',2e-1);
%             lb_stoch_prob = SReachPoint('term','voronoi-affine', sys,...
%                 initial_state, safety_tube, options);
%         end 
    end
    
    methods (Static)
        function [sys, safety_tube, initial_state] = getDI()
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
                    5e-3*eye(2)));

            % target_tube = {K, K, K, K, K, K};
            safety_tube = Tube('viability', safe_set, time_horizon);
            
            % Initial state
            initial_state = [0.4;0.4];
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
