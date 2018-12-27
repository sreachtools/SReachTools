classdef computeNormCdfInvOverApproxTest < matlab.unittest.TestCase

    methods (Test)
        function testComputeNormCdfInvOverApprox(testCase)
            n_lin_consts = 10;
            max_delta = 0.2;
    
            for pwa_accuracy = [1e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2]
                [invcdf_approx_m, invcdf_approx_c, lb_delta, norminv_knots] =...
                    computeNormCdfInvOverApprox(max_delta, pwa_accuracy, ...
                        n_lin_consts);

                x = 0;
                x(1) = lb_delta;
                for x_val = norminv_knots(2:end)
                    x = [x, linspace(x(end), x_val, 100)];
                end
                y_pwa = max(invcdf_approx_m*x+invcdf_approx_c);
                y_true = norminv(1-x);
                err_bn_pwa_and_true = y_pwa-y_true;
%                 %% Skipping plotting commands
%                 figure(); 
%                 clf
%                 plot(x,err_bn_pwa_and_true,'ro-'); hold on;
%                 plot(x,pwa_accuracy*ones(length(x),1)); 
%                 xlim([-0.05,max_delta + 0.05]);
%                 xlabel('x'); 
%                 ylabel('PWL-true curve');
%                 title('Approximation quality');
                testCase.verifyLessThanOrEqual(max(err_bn_pwa_and_true), ...
                    pwa_accuracy, 'Overapproximation error not satisfied');
                testCase.verifyLessThanOrEqual(-1e-8, ...
                    min(err_bn_pwa_and_true), ...
                    'It is not an overapproximation error');
            end
        end
    end
end
