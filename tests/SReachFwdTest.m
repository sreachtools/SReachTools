classdef SReachFwdTest < matlab.unittest.TestCase
% Unit tests for SReachFwd
% ===========================================================================
%
% This function is part of the Stochastic Optimal Control Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
%
%

    methods (Test)
        function testInputHandling(testCase)
            init_state = zeros(2,1);
            init_state_rand = zeros(2,1);
            target_time = 5;
            %% Check if state-stoch has proper input handling
            prob_str_test = 'state-stoch';
            % sys has to be a Gaussian system
            sys_uncertain = getChainOfIntegLtiSystem(2, 0.1, Polyhedron(),...
                Polyhedron('lb',-ones(2,1),'ub',ones(2,1)), eye(2));
            testCase.verifyError(@() SReachFwd(prob_str_test, sys_uncertain, init_state, target_time),'SReachTools:invalidArgs');            
            testCase.verifyError(@() SReachFwd(prob_str_test, sys_uncertain, init_state_rand, target_time),'SReachTools:invalidArgs');            
            % sys has to be uncontrolled
            sys_control = getChainOfIntegLtiSystem(2, 0.1, Polyhedron('lb',-1,'ub',1), RandomVector('Gaussian',zeros(2,1),eye(2)));
            testCase.verifyError(@() SReachFwd(prob_str_test, sys_control, zeros(2,1), target_time),'SReachTools:invalidArgs');            
            % sys has to be uncontrolled Gaussian LtiSystem
            sys = getChainOfIntegLtiSystem(2, 0.1, Polyhedron(), RandomVector('Gaussian',zeros(2,1),eye(2)));
            % Valid prob_str
            testCase.verifyError(@() SReachFwd('statestoch', sys, init_state, target_time),'SReachTools:invalidArgs');            
            % Initial state is false
            testCase.verifyError(@() SReachFwd(prob_str_test, sys, zeros(3,1), target_time),'SReachTools:invalidArgs');            
            testCase.verifyError(@() SReachFwd(prob_str_test, sys, Polyhedron(), target_time),'SReachTools:invalidArgs');            
            testCase.verifyError(@() SReachFwd(prob_str_test, sys, RandomVector('Gaussian',zeros(3,1),eye(3)), target_time),'SReachTools:invalidArgs');            
            % target_time is false
            testCase.verifyError(@() SReachFwd(prob_str_test, sys, init_state, 0),'SReachTools:invalidArgs');            
            testCase.verifyError(@() SReachFwd(prob_str_test, sys, init_state, -10),'SReachTools:invalidArgs');            
            testCase.verifyError(@() SReachFwd(prob_str_test, sys, init_state, Polyhedron()),'SReachTools:invalidArgs');            
            testCase.verifyError(@() SReachFwd(prob_str_test, sys, init_state_rand, 0),'SReachTools:invalidArgs');            
            testCase.verifyError(@() SReachFwd(prob_str_test, sys, init_state_rand, -10),'SReachTools:invalidArgs');            
            testCase.verifyError(@() SReachFwd(prob_str_test, sys, init_state_rand, Polyhedron()),'SReachTools:invalidArgs');            
            %% Probability: Assumes above input handling works
            prob_str_test = 'state-prob';
            % Empty target set
            testCase.verifyError(@() SReachFwd(prob_str_test, sys, init_state, target_time, Polyhedron.emptySet(sys.state_dim),1e-3),'SReachTools:invalidArgs');
            % Empty target tube
            testCase.verifyError(@() SReachFwd(prob_str_test, sys, init_state, target_time, Tube(Polyhedron()),1e-3),'SReachTools:invalidArgs');
            % Incorrect dimension target set
            testCase.verifyError(@() SReachFwd(prob_str_test, sys, init_state, target_time, Polyhedron('lb',-2*ones(3,1),'ub',ones(3,1)), 1e-3),'SReachTools:invalidArgs');            
            % Too strict desired_accuracy
            %testCase.verifyWarning(@() SReachFwd(prob_str_test, sys, init_state, target_time, Polyhedron('lb',-2*ones(2,1),'ub',ones(2,1)), 1e-10),'SReachTools:desiredAccuracy');            
            % More arguments than necessary
            testCase.verifyError(@() SReachFwd(prob_str_test, sys, init_state, target_time, Polyhedron('lb',-2*ones(2,1),'ub',ones(2,1)), 1e-10,'w'),'SReachTools:invalidArgs');            
             %% Probability: Assumes above input handling works
            prob_str_test = 'concat-prob';
            % Empty target set
            testCase.verifyError(@() SReachFwd(prob_str_test, sys, init_state, target_time, Polyhedron.emptySet(sys.state_dim),1e-3),'SReachTools:invalidArgs');
            % Empty target tube
            testCase.verifyError(@() SReachFwd(prob_str_test, sys, init_state, target_time, Tube(Polyhedron()),1e-3),'SReachTools:invalidArgs');
            % Incorrect dimension target set
            testCase.verifyError(@() SReachFwd(prob_str_test, sys, init_state, target_time, Tube('viability',Polyhedron('lb',-2*ones(3,1),'ub',ones(3,1)),target_time), 1e-3),'SReachTools:invalidArgs');            
            % Target tube not long enough
            testCase.verifyError(@() SReachFwd(prob_str_test, sys, init_state, target_time, Tube('viability',Polyhedron('lb',-2*ones(2,1),'ub',ones(2,1)),target_time-1), 1e-3),'SReachTools:invalidArgs');            
            % Target tube long enough --- blind test
            SReachFwd(prob_str_test, sys, init_state, target_time, Tube('viability',Polyhedron('lb',-2*ones(2,1),'ub',ones(2,1)),target_time+1), 1e-3);            
            % Too strict desired_accuracy TODO --- takes a lot of time
            %testCase.verifyWarning(@() SReachFwd(prob_str_test, sys, init_state, target_time, Polyhedron('lb',-2*ones(2,1),'ub',ones(2,1)), 1e-10),'SReachTools:desiredAccuracy');            
            % More arguments than necessary
            testCase.verifyError(@() SReachFwd(prob_str_test, sys, init_state, target_time, Tube('viability',Polyhedron('lb',-2*ones(2,1),'ub',ones(2,1)),target_time), 1e-10,'w'),'SReachTools:invalidArgs');                        
        end
        function testTerminalStateStochAndProb(testCase)
            % Performs forward stochastic reachability on the spacecraft rendezvous problem.
            % Computes the mean, covariance, and probability of lying in a
            % terminal set and compares it with Monte-Carlo simulation
            n_mcarlo_sims = 1e5;
            %% Problem setup
            % System definition
            umax=Inf;
            mean_disturbance = zeros(4,1);
            covariance_disturbance = diag([1e-4, 1e-4, 5e-8, 5e-8]);
            % Define the CWH (planar) dynamics of the deputy spacecraft relative to the chief spacecraft as a LtiSystem object
            sys_CWH = getCwhLtiSystem(4,...
                                      Polyhedron('lb', -umax*ones(2,1),...
                                                 'ub',  umax*ones(2,1)),...
                                      RandomVector('Gaussian',...
                                                   mean_disturbance,...
                                                   covariance_disturbance));
            % Create a discrete-time LQR controller that regulates the deputy to the origin 
            K = lqr(ss(sys_CWH.state_mat,sys_CWH.input_mat,[],[],-1),0.01*eye(4),eye(2));
            % Reuse the system definition in sys_CWH with appropriately defined state matrix
            closed_loop_state_mat = sys_CWH.state_mat - sys_CWH.input_mat*K;
            sys = LtiSystem('StateMatrix', closed_loop_state_mat,...
                            'DisturbanceMatrix', sys_CWH.dist_mat,...
                            'Disturbance', sys_CWH.dist);
            
            % Problem definition
            target_time = 23;                                        % Time of interest
            target_set = Polyhedron('lb',-0.05 * ones(4,1),...
                                    'ub', 0.05 * ones(4,1));          % Target set definition
            desired_accuracy = 1e-2;
            
            %% Deterministic initial state
            init_state = [-10;10;0;0]; 
            % Testing
            % Compute mean and cov
            [mean_vec, cov_mat] = SReachFwd('state-stoch', sys, init_state, target_time);
            % Compute probability
            prob = SReachFwd('state-prob', sys, init_state, target_time, target_set, desired_accuracy);
            % Validation using Monte-Carlo
            % This function returns the concatenated state vector stacked columnwise
            concat_state_realization = generateMonteCarloSims(...
                                                           n_mcarlo_sims,...
                                                           sys,...
                                                           init_state,...
                                                           target_time);
            % Extract the location of the deputy at target_time
            end_locations = concat_state_realization(end-sys.state_dim +1 : end,:);
            % Check if the location is within the target_set or not
            mcarlo_result = target_set.contains(end_locations);
            mean_end_locations = mean(end_locations,2);
            % Checking values
            testCase.verifyLessThanOrEqual(abs(sum(mcarlo_result)/n_mcarlo_sims - prob), desired_accuracy);            
            testCase.verifyLessThanOrEqual(max(abs(mean_vec - mean_end_locations)), 1e-2);            
            %% Stochastic initial state
            init_state_rand = RandomVector('Gaussian',init_state,0.1*eye(4)); 
            % Testing
            % Compute mean and cov
            [mean_vec, cov_mat] = SReachFwd('state-stoch', sys, init_state_rand, target_time);
            % Compute probability
            prob = SReachFwd('state-prob', sys, init_state_rand, target_time, target_set, desired_accuracy);
            % Validation using Monte-Carlo
            % This function returns the concatenated state vector stacked columnwise
            concat_state_realization = generateMonteCarloSims(...
                                                           n_mcarlo_sims,...
                                                           sys,...
                                                           init_state_rand,...
                                                           target_time);
            % Extract the location of the deputy at target_time
            end_locations = concat_state_realization(end-sys.state_dim +1 : end,:);
            % Check if the location is within the target_set or not
            mcarlo_result = target_set.contains(end_locations);
            mean_end_locations = mean(end_locations,2);
            % Checking values
            testCase.verifyLessThanOrEqual(abs(sum(mcarlo_result)/n_mcarlo_sims - prob), desired_accuracy);            
            testCase.verifyLessThanOrEqual(max(abs(mean_vec - mean_end_locations)), 1e-2);                    
        end
        
        function testTubeStochAndProb(testCase)
            % Performs forward stochastic reachability on the spacecraft rendezvous problem.
            % Computes the mean, covariance, and probability of lying in a
            % terminal set and compares it with Monte-Carlo simulation
            n_mcarlo_sims = 1e5;
            %% Problem setup
            % System definition
            umax=Inf;
            mean_disturbance = zeros(4,1);
            covariance_disturbance = diag([1e-4, 1e-4, 5e-8, 5e-8]);
            % Define the CWH (planar) dynamics of the deputy spacecraft relative to the chief spacecraft as a LtiSystem object
            sys_CWH = getCwhLtiSystem(4,...
                                      Polyhedron('lb', -umax*ones(2,1),...
                                                 'ub',  umax*ones(2,1)),...
                                      RandomVector('Gaussian',...
                                                   mean_disturbance,...
                                                   covariance_disturbance));
            % Create a discrete-time LQR controller that regulates the deputy to the origin 
            K = lqr(ss(sys_CWH.state_mat,sys_CWH.input_mat,[],[],-1),0.01*eye(4),eye(2));
            % Reuse the system definition in sys_CWH with appropriately defined state matrix
            closed_loop_state_mat = sys_CWH.state_mat - sys_CWH.input_mat*K;
            sys = LtiSystem('StateMatrix', closed_loop_state_mat,...
                            'DisturbanceMatrix', sys_CWH.dist_mat,...
                            'Disturbance', sys_CWH.dist);
            
            % Problem definition
            target_time = 10;                                        % Time of interest
            target_set = Polyhedron('lb',-0.05 * ones(4,1),...
                                    'ub', 0.05 * ones(4,1));         % Target set definition
            %% Safe set definition --- LoS cone |x|<=y and y\in[0,ymax] and |vx|<=vxmax and |vy|<=vymax
            ymax=10;
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
            % Create a target tube
            target_tube = Tube('reach-avoid', safe_set, target_set, target_time);
            desired_accuracy = 1e-2;
            
            %% Deterministic initial state
            init_state = [0;-1;0;0]; 
            % Testing
            % Compute mean and cov
            [mean_vec, cov_mat] = SReachFwd('concat-stoch', sys, init_state, target_time);
            % Compute probability
            prob = SReachFwd('concat-prob', sys, init_state, target_time, target_tube, desired_accuracy);
            % Validation using Monte-Carlo
            % This function returns the concatenated state vector stacked columnwise
            concat_state_realization = generateMonteCarloSims(...
                                                           n_mcarlo_sims,...
                                                           sys,...
                                                           init_state,...
                                                           target_time);
            % Check if the location is within the target_set or not
            mcarlo_result = target_tube.contains([repmat(init_state,1,n_mcarlo_sims);
                                      concat_state_realization]);
            mean_concat_state = mean(concat_state_realization,2);
            % Checking values
            testCase.verifyLessThanOrEqual(abs(sum(mcarlo_result)/n_mcarlo_sims - prob), desired_accuracy);            
            testCase.verifyLessThanOrEqual(max(abs(mean_vec - mean_concat_state)), 1e-2);            
            %% Stochastic initial state
            init_state_rand = RandomVector('Gaussian',init_state,0.001*eye(4)); 
            % Testing
            % Compute mean and cov
            [mean_vec, cov_mat] = SReachFwd('concat-stoch', sys, init_state_rand, target_time);
            % Compute probability
            prob = SReachFwd('concat-prob', sys, init_state_rand, target_time, target_tube, desired_accuracy);
            % Validation using Monte-Carlo
            % This function returns the concatenated state vector stacked columnwise
            concat_state_realization = generateMonteCarloSims(...
                                                           n_mcarlo_sims,...
                                                           sys,...
                                                           init_state_rand,...
                                                           target_time);
            % Check if the location is within the target_set or not
            mcarlo_result = target_tube.contains([repmat(init_state,1,n_mcarlo_sims);
                                      concat_state_realization]);
            mean_concat_state = mean(concat_state_realization,2);
            % Checking values
            testCase.verifyLessThanOrEqual(abs(sum(mcarlo_result)/n_mcarlo_sims - prob), desired_accuracy);            
            testCase.verifyLessThanOrEqual(max(abs(mean_vec - mean_concat_state)), 1e-2);                    
        end
    end
end

