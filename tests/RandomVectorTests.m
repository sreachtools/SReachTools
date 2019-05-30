classdef RandomVectorTests < matlab.unittest.TestCase
% SReachTools/TubeTests: Unit tests for bounded disturbances
% ===========================================================================
%
% Unit testing for RandomVector
%
% ===========================================================================
%
% This function is part of the Stochastic Optimal Control Toolbox.
% License for the use of this function is given in
%      https://sreachtools.github.io/license/
%
%

    methods (Test)
        function inputHandlingTest(test_case)
            % Gaussian: Too many arguments => Varargin means invalidArgs
            test_case.verifyError(@(x) RandomVector('Gaussian',zeros(2,1),...
                eye(2), eye(2)), 'SReachTools:invalidArgs');
            % Gaussian: Not enough arguments
            test_case.verifyError(@(x) RandomVector('Gaussian',zeros(2,1)),...
                'SReachTools:invalidArgs');
            % Gaussian: Input parsing
            test_case.verifyError(@(x) RandomVector('Gaussian',zeros(2,1),...
                ones(2,3)), 'MATLAB:RandomVector:RandomVector:expectedSquare');
            test_case.verifyError(@(x) RandomVector('Gaussian',zeros(2,1),...
                ones(4,3)), 'MATLAB:RandomVector:RandomVector:expectedSquare');
            test_case.verifyError(@(x) RandomVector('Gaussian',zeros(1,3),...
                eye(3)), 'MATLAB:RandomVector:RandomVector:expectedColumn');
            % Gaussian: Dimension mismatch
            test_case.verifyError(@(x) RandomVector('Gaussian',zeros(2,1),...
                eye(3)), 'SReachTools:invalidArgs');
            % Gaussian: Non-symmetric matrix
            warning('on','SReachTools:runtime');
            test_case.verifyWarning(@(x) RandomVector('Gaussian',zeros(2,1),...
                eye(2)+[0,3e-16;-3e-16,0]), 'SReachTools:runtime');
            RandomVector('Gaussian',zeros(2,1), eye(2)+[0,1e-18;-1e-18,0]);
            % Gaussian: Non-positive semi-definite matrix
            test_case.verifyError(@(x) RandomVector('Gaussian',zeros(2,1),...
                [-1,0;0,0]), 'SReachTools:invalidArgs');
            %% Gaussian: Deterministic Gaussian semi-definite matrix
            %test_case.verifyWarning(@(x) RandomVector('Gaussian',zeros(2,1),...
            %    [1,0;0,0]), 'SReachTools:runtime');
            % Invalid string
            test_case.verifyError(@(x) RandomVector('Exp',zeros(2,1)),...
                'MATLAB:unrecognizedStringChoice');
            % Simply define a well-defined Gaussian
            RandomVector('Gaussian',zeros(2,1), eye(2));
            % Simply define a well-defined UserDefined RV: exp
            RandomVector('UserDefined',@(N) exprnd(1,N));
        end
        function multiplicationTest(test_case)        
            %% Define a well-defined Gaussian
            r = RandomVector('Gaussian',zeros(2,1), eye(2));
            eye(2) * r;
            r * eye(2);
            [1, 0 ] * r;
            %warning('on','SReachTools:runtime');
            %test_case.verifyWarning(@(x) [1, 0;1, 0;1, 0] * r,...
            %    'SReachTools:runtime');
            % Invalid dimension
            test_case.verifyError(@(x) r * eye(3), 'SReachTools:invalidArgs');            
            % Invalid input
            test_case.verifyError(@(x) r * 'ch', 'SReachTools:invalidArgs');            

            % Repeat this for exponential case (generator)
            exp_mean = 2;
            r_exp = RandomVector.exponential(exp_mean);
            newr_exp = 3 * r_exp;
            test_case.verifyTrue(abs(3*r_exp.mean() - newr_exp.mean()) <1e-1,...
                sprintf(['Mismatch in mean for exponential case | MC: %1.3e',...
                    ' expected: %1.3f'], newr_exp.mean(), 3 * exp_mean));
                
            % Invalid dimension
            test_case.verifyError(@(x) r_exp*eye(3), 'SReachTools:invalidArgs');            
            % Invalid input
            test_case.verifyError(@(x) r_exp * 'ch', 'SReachTools:invalidArgs');            
        end
        
        function plusTest(test_case)
            % Define a well-defined Gaussian
            r = RandomVector('Gaussian',zeros(2,1), eye(2));
            newr = ones(2,1) + r;
            test_case.verifyTrue(isequal(r.mean(), newr.mean() - ones(2,1)),...
                'Mismatch in mean');
            test_case.verifyTrue(isequal(r.cov(), newr.cov()),...
                'Mismatch in covariance');
            % Invalid dimension
            test_case.verifyError(@(x) r + ones(3,1),'SReachTools:invalidArgs');            
            % Invalid input
            test_case.verifyError(@(x) r + 'ch', 'SReachTools:invalidArgs');            

            %% Repeat this for exponential case (generator)
            r_exp = RandomVector.exponential(2);
            newr_exp = 3 + r_exp;
            test_case.verifyTrue(abs(3+r_exp.mean() - newr_exp.mean())<1e-2,...
                'Mismatch in mean for exponential case');
            test_case.verifyError(@(x) [1, 0;1, 0;1, 0] + r_exp,...
                'SReachTools:invalidArgs');            
            % Invalid dimension
            test_case.verifyError(@(x) r_exp + ones(3,1),...
                'SReachTools:invalidArgs');            
            % Invalid input
            test_case.verifyError(@(x) r_exp + 'ch', 'SReachTools:invalidArgs');            
            
            %% Gaussian + F*Exp
            gauss_plus_exp = ones(2,1) * r_exp + r;
            test_case.verifyTrue(isequal(gauss_plus_exp.type,'UserDefined'),...
                'Expected it to be UserDefined');
            test_case.verifyTrue(...
                max(abs(gauss_plus_exp.mean() -...
                    ones(2,1) * r_exp.mean()-r.mean()))<1e-2,...
                'Mismatch in mean');

            %% Exp + Exp
            exp_plus_exp = r_exp + r_exp;
            test_case.verifyTrue(isequal(exp_plus_exp.type,'UserDefined'),...
                'Expected it to be UserDefined');
            test_case.verifyTrue(...
                abs(exp_plus_exp.mean() - 2 * r_exp.mean())<1e-2,...
                'Mismatch in mean');
            
            %% Gauss + Gauss
            gauss_plus_gauss = r + r;
            test_case.verifyTrue(isequal(gauss_plus_gauss.type,'Gaussian'),...
                'Expected it to be UserDefined');
            test_case.verifyTrue(...
                max(abs(gauss_plus_gauss.mean() - 2 * r.mean()))<1e-8,...
                'Mismatch in mean');
            
        end
        function concatTest(test_case)        
            % Define a well-defined Gaussian
            r = RandomVector('Gaussian',zeros(2,1), eye(2));
            R = r.concat(10);
            test_case.verifyTrue(nnz(R.mean())==0, ['Should have ',...
                'all been zero']);
            test_case.verifyTrue(isequal(size(R.mean()), [20,1]),...
                'Dimension mismatch');
        end
        
        function meanCovStaticMethodAndgetRealizationTest(test_case)        
            % Test mean and covariance of static methods
            
            mean_vec = [3;4];
            temp_mat = rand(2,2);
            cov_matrix = [temp_mat + temp_mat']/2 + diag([1,100]);
            r = RandomVector.gaussian(mean_vec, cov_matrix);
            mean_gauss_mcarlo = mean(r.getRealizations(1e6), 2);
            test_case.verifyTrue(isequal(mean_vec, r.mean()), ['Mismatch in',...
                ' mean']);
            test_case.verifyTrue(isequal(cov_matrix, r.cov()), ['Mismatch ',...
                'in cov']);
            % Relying on Monte-Carlo: Sample mean can be off
            test_case.verifyTrue(max(abs(mean_vec - mean_gauss_mcarlo))<1e-1,...
                'Mismatch in mean (However, it is monte-carlo vs exact)');
            
            % Define an exponential distribution
            mean_vec = [3;1]; % 1/lambda
            expected_mean = mean_vec;
            expected_cov = diag(expected_mean.^2); % 1/lambda^2
            r = RandomVector.exponential(mean_vec);
            % Test r.mean() and r.cov() which in turn tests
            % getRealization() 
            % Relying on Monte-Carlo: Sample mean/cov can be off
            mc_mean = r.mean();
            test_case.verifyTrue(max(max(abs(expected_mean - mc_mean)))<1e-1,...
                sprintf('Mismatch in mean MC: %1.3e expected: %1.3f',...
                mc_mean, expected_mean));
            mc_cov = r.cov();
            test_case.verifyTrue(max(max(abs(expected_cov - mc_cov)))<1e-1,...
                sprintf('Mismatch in mean MC: %1.3e expected: %1.3f',...
                mc_cov, expected_cov));
        end
        
        function getProbPolyhedronTest(test_case)        
            % Test probability of set
            
            desired_accuracy = 1e-2;
            % Gaussian case            
            mean_vec = [3;4];
            temp_mat = rand(2,2);
            cov_matrix = [temp_mat + temp_mat']/2 + diag([1,1]);
            r = RandomVector.gaussian(mean_vec, cov_matrix);
            
            lb_polytope = mean_vec - 2*ones(2,1);
            ub_polytope = mean_vec + 2*ones(2,1);
            
            test_polyhedron = Polyhedron('lb', lb_polytope,...
                'ub', ub_polytope);
            prob = r.getProbPolyhedron(test_polyhedron, desired_accuracy);
            prob_mvncdf = mvncdf(lb_polytope', ub_polytope', mean_vec',...
                cov_matrix);
            test_case.verifyLessThanOrEqual(prob, prob_mvncdf);
            
            
            % Exponential case
            % Define an exponential distribution
            mean_vec = [1;3]; % 1/lambda
            max_corner = 2;
            r = RandomVector.exponential(mean_vec);
            lb_polytope = zeros(2,1);
            ub_polytope = max_corner *ones(2,1);
            
            test_polyhedron = Polyhedron('lb', lb_polytope,...
                'ub', ub_polytope);
            prob = r.getProbPolyhedron(test_polyhedron, desired_accuracy);            
            prob_expcdf = prod(expcdf(max_corner, mean_vec));            
            % Comparing monte carlo to true value
            test_case.verifyLessThanOrEqual(prob, prob_expcdf);
        end
        
        function vertCatTest(testCase)        
            % Test probability of set
            
            % Gaussian case            
            mean_vec = [3;4];
            temp_mat = rand(2,2);
            cov_matrix = [temp_mat + temp_mat']/2 + diag([1,1]);
            r_gauss = RandomVector.gaussian(mean_vec, cov_matrix);
            
            rnew = [r_gauss;r_gauss;r_gauss];
            testCase.verifyTrue(max(abs(rnew.mean() - ...
                [mean_vec;mean_vec;mean_vec])) ...
                <1e-8,'Failed concatenation of Gaussian random vectors (mean)');
            testCase.verifyTrue(max(max(abs(rnew.cov() - ...
                blkdiag(cov_matrix, cov_matrix, cov_matrix))))<1e-8, ...
                'Failed concatenation of Gaussian random vectors (cov)');
            testCase.verifyTrue(strcmpi(rnew.type,r_gauss.type), ...
                'Failed concatenation of Gaussian random vectors (type)');
            
            % Exponential case
            % Define an exponential distribution
            mean_vec = [1;3]; % 1/lambda
            r_exp = RandomVector.exponential(mean_vec);
            rnew = [r_exp;r_exp;r_exp];
            testCase.verifyTrue( ...
                max(abs(rnew.mean() - [mean_vec;mean_vec;mean_vec]))<1e-2, ...
                'Failed concatenation of UserDefined random vectors (mean)');
            testCase.verifyTrue(strcmpi(rnew.type,r_exp.type), ...
                'Failed concatenation of Gaussian random vectors (type)');
            
            testCase.verifyError(@() [r_gauss;r_exp],'SReachTools:invalidArgs');
        end
    end
end
