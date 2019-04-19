function [overapprox_m, overapprox_c, lb_phiinv, norminv_knots] =...
    computeNormCdfInvOverApprox(max_delta, pwa_accuracy, n_lin_consts)
% Compute a piecewise-linear overapproximation of norminv(1-x) for 
% x \in [lb_delta,0.5] to a user-specified quality (Internal function)
% =============================================================================
%
% computeNormCdfInvOverApprox generates a piecewise-linear overapproximation of
% norminv(1-x) for x\in[lb_delta,0.5]. Specifically, given any
% z\in[lb_delta,0.5],
%
% norminv(1-x) + pwa_accuracy>max(cdf_approx_m * x + cdf_approx_c)>norminv(1-x),
%
% with lb_delta = max(max_delta/n_lin_consts/10,1e-8).
%
% This function implements Algorithm 2 of the following paper:
%
% A. Vinod and M. Oishi. Affine controller synthesis for stochastic reachability
% via difference of convex programming. In Proc. Conf. Dec. & Ctrl. 2019.
% (submitted). https://hscl.unm.edu/affinecontrollersynthesis/
%
% =============================================================================
%
% 
% Inputs: 
% -------
%   max_delta    - This is the maximum tolerance for violation of the joint
%                  chance constraint | risk allocation can not exceed this value
%   pwa_accuracy - Accuracy expected from the chance constraint formulation with 
%                  n_lin_consts, no. of individual chance constraints
%   n_lin_consts - No. of individual chance constraints
%
% Outputs:
% --------
%   overapprox_m - Secant slopes that will overapproximate norminv(1-x)
%   overapprox_c - Secant y-intercepts that will overapproximate norminv(1-x)
%   lb_phiinv    - Lower bound on x in norminv(1-x) for which the provided PWA
%                   approximation
%   norminv_knots- Breakpoints of the PWA overapproximation, i.e., points
%                   at which the PWA overapproximation coincides with the 
%                   norminv curve
%
% Notes:
% * Partial INPUT HANDLING: Checks only if the bounds provided are accurate
% * MATLAB DEPENDENCY: None
% * The requested desired_accuracy/n_lin_consts/10 can not be smaller than 1e-8.
% 
% ==============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://sreachtools.github.io/license/
%
%
    
    % Interested bounds are desired_accuracy/N_{ineq} and max_delta
    lb_phiinv = max(pwa_accuracy/n_lin_consts/10,1e-8);
    ub_phiinv = max_delta;
    
    % Throw an error if the requested bounds are beyond the PWA parameters
    if ub_phiinv > 0.5
        throwAsCaller(SrtInvalidArgsError(['Upper bound can''t exceed 0.5 ', ...
            '(breaks convexity).']));
    end
    
    if pwa_accuracy < 1e-5
        warning('SReachTools:desiredAccuracy', ...
            ['The requested accuracy might take a lot of time, and may ', ...
             ' cause MATLAB to crash.']);
    end
    
    %% hessian_phinvOneMinusX definition: Hessian of phiinv(1-x)
    % See end of the file for MATLAB's symbolic toolbox commands to obtain
    % this function
    hessian_phiinvOneMinusX = @(x) 2*2^(1/2)*pi*erfcinv(2*x) * ...
        exp(2*erfcinv(2*x)^2);
    
	%% Initialization
    norminv_knots(1) = lb_phiinv;
    j = 1;
    
    % Iterate till we reach the end
    while norminv_knots(j)< ub_phiinv
        hval = sqrt(8*pwa_accuracy/hessian_phiinvOneMinusX(norminv_knots(j)));
        norminv_knots(j+1) = min([norminv_knots(j) + hval,ub_phiinv]);
        % disp(knots_underapprox(j));
        
        %% Construction of the underapproximation via Lagrange interpolation
        y_2 = norminv(1-norminv_knots(j+1));
        y_1 = norminv(1-norminv_knots(j));
        x_2 = norminv_knots(j+1);
        x_1 = norminv_knots(j);

        % Lagrange linear interpolation
        overapprox_m(j,1) = (y_2 - y_1)/(x_2 - x_1);
        overapprox_c(j,1) = y_1 - overapprox_m(j) * x_1;
        
        %% Increment j
        j = j+1;       
    end
end
