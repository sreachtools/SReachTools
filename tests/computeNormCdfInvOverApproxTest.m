classdef computeNormCdfInvOverApproxTest < matlab.unittest.TestCase

    methods (Test)
        function testComputeNormCdfInvOverApproxExceedsLimits(testCase)

            desired_accuracy = 1e-6;
            n_lin_consts = 200;
            max_delta = 0.2;

            testCase.verifyError(@() computeNormCdfInvOverApprox(max_delta, desired_accuracy,...
                    n_lin_consts),'SReachTools:invalidArgs');
        end

        function testComputeNormCdfInvOverApprox(testCase)

            desired_accuracy = 1e-4;
            n_lin_consts = 100;
            max_delta = 0.2;

            [cdf_approx_m, cdf_approx_c, lb_delta] =...
                computeNormCdfInvOverApprox(max_delta, desired_accuracy,...
                    n_lin_consts);
                
            x_pwl = [lb_delta:1e-8:1e-6,...
                     1e-6:1e-7:1e-5,...
                     1e-5:1e-6:1e-4,...
                     1e-4:1e-5:1e-3,...
                     1e-3:1e-4:max_delta];
            x(1) = lb_delta;
            for x_val = x_pwl(2:end)
                x = [x, linspace(x(end), x_val, 100)];
            end
            y=max(cdf_approx_m*x+cdf_approx_c);
            y_true = norminv(1-x);
            err=y-y_true;
            
            testCase.verifyLessThanOrEqual(lb_delta, desired_accuracy/n_lin_consts/10,'Incorrect lower bound')
            testCase.verifyLessThanOrEqual(max(err), desired_accuracy,'Larger than allowed error')

            %% Skipping plotting commands
%             figure(); plot(x,y,'ro-'); hold on; plot(x,y_true,'b*-'); legend('PWL','True')
%             xlabel('x'),ylabel('norminv(1-x)');title('Approximation');
%             figure(); plot(x,err,'ro-');hold on;plot(x,desired_accuracy*ones(length(x),1)); 
%             xlim([-0.05,max_delta + 0.05]);ylim([-1e-5,1.5*desired_accuracy]);  xlabel('x'); ylabel('PWL-true curve');
%             title('Approximation quality');
        end
        
        function testComputeNormCdfInvOverApproxAfterCreationExceedsLimits(testCase)

            desired_accuracy = 1e-8;
            n_lin_consts = 200;
            max_delta = 0.2;

            testCase.verifyError(@() computeNormCdfInvOverApprox(max_delta, desired_accuracy,...
                    n_lin_consts),'SReachTools:invalidArgs');
        end
    end
end
