function [cdf_approx_A, cdf_approx_b, lower_bound] =...
    computeNormCdfInvOverApprox(no_samples)
% Computes the piece-linear overapproximation of norm-inv(1-x) from [0,0.5]

    x = linspace(0,0.5,no_samples);
    x = x(2:end);
    lower_bound = x(2);
    
    cdf_approx_A = zeros(no_samples - 2, 1);
    cdf_approx_b = zeros(no_samples - 2, 1);
    
    for indx_x = 1:no_samples-2
        y_2 = norminv(1-x(indx_x+1));
        y_1 = norminv(1-x(indx_x));
        x_2 = x(indx_x + 1);
        x_1 = x(indx_x);
        cdf_approx_A(indx_x) = (y_2 - y_1)/(x_2 - x_1);
        cdf_approx_b(indx_x) = y_1 - cdf_approx_A(indx_x) * x_1;
    end    
end
