function [cdf_approx_m, cdf_approx_c,varargout] = computeNormCdfInvOverApprox(varargin)
% Compute a piecewise-linear overapproximation of norminv(1-x) for 
% x \in [1e-5,0.5] to the quality of 1e-4
% =============================================================================
%
% computeNormCdfInvOverApprox generates a piecewise-linear overapproximation of
% norminv(1-x) for x\in[1e-5,0.5]. Specifically, given any z\in[1e-5,0.5],
% norm(1-z) +1e-4 > max(cdf_approx_m * z + cdf_approx_c) > norm(1-z).
%
% USAGE: See getUnderapproxStochReachAvoidSet,
% computeCcLowerBoundStochReachAvoidPwlRisk.
%
% =============================================================================
%
% 
% Inputs: None
% -------
%
% Outputs:
% --------
%   cdf_approx_m - Secant slopes that will overapproximate norminv(1-x)
%   cdf_approx_c - Secant y-intercepts that will overapproximate norminv(1-x)
%   lb_x         - (Optional) Lower-bound on the range of the approximation
%   diff_val     - (Optional) Smallest step-size between the secant end points
%
% Notes:
% * NOT ACTIVELY TESTED: TODO
% * NO INPUT HANDLING: No arguments needed
% * MATLAB DEPENDENCY: None
% * Max error of approximation is 0.99818e-04 (estimated to a nbd of 1e-7)
% * The end points of the secants are obtained by the sequence 
%       {lb_x, lb_x + h gamma^{0:n_x},0.5, and the midpoints} 
%   where lb_x =1e-5, h=1e-6, gamma=1.088, n_x =|_log((0.5-lb_x)/h)/log(gamma)_|
% * If we desire to optimize delta_i where i\in [1,M], then this approach
%   results in a maximum artificial conservativeness (on top of Boolean
%   conservativeness) of 1e-4*M.
% 
% ================= ============================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
%
%
    
    % Compute the sequence {lb_x, lb_x + h *gamma^{0:n_x},0.5}
    lb_x = 1e-5;
    if nargin == 0 || abs(varargin{1}-1e-4)<eps
        max_error_estim = 1e-4;
        h = 1e-6;
        gamma = 1.088;
    else
        max_error_estim = 1e-3;
        h = 2e-6;
        gamma = 1.29;
    end
    n_x = floor(log((0.5 - lb_x)/h)/log(gamma));
    x = lb_x + h.*(gamma.^(0:n_x));
    x = [lb_x,x,0.5];

    % Compute the midpoints of the sequence
    x_1 = x(1:end-1)/2 + x(2:end)/2;    

    % Add to the existing list and sort it
    x = sort([x,x_1]);                              %gamma = 1.088; n=314

    % no. of lines joining n points is n-1
    n_piecewise_lin_comp = length(x)-1;

    % Initialize the vectors for secant information
    cdf_approx_m = zeros(n_piecewise_lin_comp, 1);
    cdf_approx_c = zeros(n_piecewise_lin_comp, 1);

    % Compute the secants    
    for indx_x = 1:n_piecewise_lin_comp
        y_2 = norminv(1-x(indx_x+1));
        y_1 = norminv(1-x(indx_x));
        x_2 = x(indx_x + 1);
        x_1 = x(indx_x);
        % y=mx+c where m and c are computed from secant end points
        cdf_approx_m(indx_x) = (y_2 - y_1)/(x_2 - x_1);
        cdf_approx_c(indx_x) = y_1 - cdf_approx_m(indx_x) * x_1;
    end    
    varargout{1} = lb_x;
    varargout{2} = max_error_estim;     % Max error estimate
    varargout{3} = x;
end

%% Other options
%     x_1 = x(1:end-1)/4 + x(2:end)*3/4;
%     x_2 = x(1:end-1)/2 + x(2:end)/2;
%     x_3 = x(1:end-1)*3/4 + x(2:end)/4; 
%     x = sort([x,x_1,x_2,x_3]);                      %gamma = 1.175; n=332
%     x_1 = x(1:end-1)/3 + x(2:end)*2/3;
%     x_2 = x(1:end-1)*2/3 + x(2:end)/3;
%     x = sort([x,x_1,x_2]);                          %gamma = 1.13;  n=327


%% Sanity check
% See the quality of the piecewise linear approximation in the range [1e-5,0.5]
% by running the following command
%
% clear;close all;
% [cdf_approx_m, cdf_approx_c,a, max_error_pwl, x_pwl] = computeNormCdfInvOverApprox();
% x = x_pwl(1);
% for x_val = x_pwl(2:end)
%     x = [x, linspace(x(end), x_val, 100)];
% end
% y=max(cdf_approx_m*x+cdf_approx_c);
% y_true = norminv(1-x);
% err=y-y_true;
% figure(); plot(x,y,'ro-'); hold on; plot(x,y_true,'b*-'); legend('PWL','True')
% xlabel('x'),ylabel('norminv(1-x)');title('Approximation');
% figure(); plot(x,err,'ro-');hold on;plot(x,1e-4*ones(length(x),1)); 
% xlim([-0.2,0.5]);ylim([-1e-5,1.5e-4]);  xlabel('x'); ylabel('PWL-true curve');
% title('Approximation quality');
% if max(err) > max_error_pwl
%   disp('Predicted error is smaller than the original error');
% end
% fprintf(['Min error: %1.4e | Max error: %1.4e | No. of ineq: %d\n'],...
%   min(err), max(err), length(cdf_approx_m));
