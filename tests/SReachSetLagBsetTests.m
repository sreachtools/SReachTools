classdef SReachSetLagBsetTests < matlab.unittest.TestCase
% Socbox/BoundedDisturbanceTests: Unit tests for bounded disturbances
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
%      https://github.com/unm-hscl/Socbox/blob/master/LICENSE
%
%

    methods (Test)
        
        function testInputHandling(test_case)
        end
        
        function testBsetGeneration(test_case)
        % Socbox/BoundedDisturbanceTests/testBoundedEllipseRandomDirections: 
        % Unit test for getBsetForDisturbance for random ellipsoid 
        % generation
        % ====================================================================
        %
        % Unit test for getBsetForDisturbance using ellipsoid 
        % overapproximation via random direction selection
        %   
        % ====================================================================
        %
        % This function is part of the Stochastic Optimal Control Toolbox.
        % License for the use of this function is given in
        %      https://github.com/unm-hscl/Socbox/blob/master/LICENSE
        %
        %   

            T = 0.25;
            umax = 0.75;

            sys = LtiSystem(...
                'StateMatrix', [1, T; 0, 1], ...
                'InputMatrix', [T^2; T], ...
                'InputSpace', Polyhedron('lb', -umax, 'ub', umax), ...
                'DisturbanceMatrix', eye(2), ...
                'Disturbance', RandomVector('Gaussian',zeros(2,1),0.01*eye(2)));

            prob_thresh = 0.7;
            time_horizon = 5;
    
            %% Random sampling of the ellipse
            options = SReachSetOptions('term','lag-under', ...
                'bound_set_method','random', ...
                'num_dir',100);
            bounded_dist = SReachSetLagBset(sys, ...
                prob_thresh^(1/time_horizon), options);
            test_case.verifyLessThanOrEqual(abs(prob_thresh^(1/time_horizon)...
                - test_case.computeProb(sys.dist, bounded_dist)), 1e-2);
            test_case.verifyInstanceOf(bounded_dist, 'Polyhedron');
            
%             %% Optim-box generation
%             options = SReachSetOptions('term','lag-under', ...
%                 'bound_set_method','optim-box', ...
%                 'box_center',ones(2,1));
%             bounded_dist = SReachSetLagBset(sys.dist, time_horizon, ...
%                 prob_thresh, options);      
%             test_case.verifyLessThanOrEqual(abs(prob_thresh^(1/time_horizon)...
%                 - test_case.computeProb(sys.dist, bounded_dist)), 1e-2);            
%             test_case.verifyInstanceOf(bounded_dist, 'Polyhedron');
            
%             %% Loading option
%             options = SReachSetOptions('term','lag-under', ...
%                 'bound_set_method','load', ...
%                 'load_str','data/cwhUnderapproxBoundeSet.mat');
%             bounded_dist = SReachSetLagBset([],[],[], options);            
%             
%             dist = RandomVector('Gaussian', zeros(4,1), ...
%                 diag([1e-4, 1e-4, 5e-8, 5e-8]));
%             test_case.verifyInstanceOf(bounded_dist, 'Polyhedron');           
%             test_case.verifyLessThanOrEqual(abs(prob_thresh^(1/time_horizon)...
%                 - test_case.computeProb(dist, bounded_dist)), 1e-2);
            
%             %% Box generation
%             options = SReachSetOptions('term','lag-under', ...
%                 'bound_set_method','box', ...
%                 'err_thresh',1e-2);
%             bounded_dist = SReachSetLagBset(sys.dist, time_horizon, ...
%                 prob_thresh, options);            
%             test_case.verifyInstanceOf(bounded_dist, 'Polyhedron');           
%             test_case.verifyLessThanOrEqual(abs(prob_thresh^(1/time_horizon)...
%                 - test_case.computeProb(sys.dist, bounded_dist)), 1e-2);
        end
    end
    
    methods (Static)
        function prob = computeProb(dist, bset)
            % Construct the half-space representation for qscmvnv
            cov_mat = (dist.parameters.covariance + ...
                dist.parameters.covariance')/2; 
            qscmvnv_lb = repmat(-Inf, [size(bset.A, 1), 1]);
            qscmvnv_coeff_matrix = bset.A;
            qscmvnv_ub = bset.b - bset.A * dist.parameters.mean;
            prob = iteratedQscmvnv(cov_mat, ...
                                   qscmvnv_lb, ...
                                   qscmvnv_coeff_matrix, ...
                                   qscmvnv_ub, ...
                                   1e-3, ...
                                   10);            
        end
    end
end
