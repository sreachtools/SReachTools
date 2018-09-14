classdef computeNormCdfInvOverApproxTest < matlab.unittest.TestCase

    methods (Test)
        function testComputeNormCdfInvOverApprox(testCase)
            % Warning is stored, disabled and enabled based on need, and
            % then restored
            n_lin_consts = 10;
            max_delta = 0.2;
            err_ranges = [1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2];
    
            orig_warning_state = warning;
            for desired_accuracy = [1e-6, 1e-5, 2e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2]
                warning('off','all'); 
                [cdf_approx_m, cdf_approx_c, lb_delta, norminv_knots] =...
                    computeNormCdfInvOverApprox(max_delta, desired_accuracy,...
                        n_lin_consts);
                warning('on','all');
                
                x(1) = lb_delta;
                for x_val = norminv_knots(2:end)
                    x = [x, linspace(x(end), x_val, 100)];
                end
                y=max(cdf_approx_m*x+cdf_approx_c);
                y_true = norminv(1-x);
                err=y-y_true;
                
                testCase.verifyLessThanOrEqual(lb_delta, desired_accuracy/n_lin_consts/10,'Incorrect lower bound')
                if desired_accuracy < err_ranges(1)
                    testCase.verifyWarning(@() computeNormCdfInvOverApprox(max_delta, desired_accuracy,...
                        n_lin_consts),'','');                
                    testCase.verifyLessThanOrEqual(max(err), err_ranges(1),'Larger than allowed error');
                else
                    testCase.verifyLessThanOrEqual(max(err), desired_accuracy,'Larger than allowed error');
                end
            end
            % Too far away
            warning('off','all');                 
            testCase.verifyError(@() computeNormCdfInvOverApprox(max_delta, 1e-7,...
                        n_lin_consts),'SReachTools:invalidArgs');
            warning(orig_warning_state);                 
            %% Skipping plotting commands
%             figure(); plot(x,y,'ro-'); hold on; plot(x,y_true,'b*-'); legend('PWL','True')
%             xlabel('x'),ylabel('norminv(1-x)');title('Approximation');
%             figure(); plot(x,err,'ro-');hold on;plot(x,desired_accuracy*ones(length(x),1)); 
%             xlim([-0.05,max_delta + 0.05]);ylim([-1e-5,1.5*desired_accuracy]);  xlabel('x'); ylabel('PWL-true curve');
%             title('Approximation quality');
        end
    end
end
