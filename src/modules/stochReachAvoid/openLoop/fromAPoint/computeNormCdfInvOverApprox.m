function [cdf_approx_m, cdf_approx_c, lb_phiinv] =...
    computeNormCdfInvOverApprox(max_delta, desired_accuracy, n_lin_consts)
% Compute a piecewise-linear overapproximation of norminv(1-x) for 
% x \in [1e-5,0.5] to the quality of 1e-4
% =============================================================================
%
% computeNormCdfInvOverApprox generates a piecewise-linear overapproximation of
% norminv(1-x) for x\in[1e-5,0.5]. Specifically, given any z\in[1e-5,0.5],
% norminv(1-x) + err_bnd > max(cdf_approx_m * x + cdf_approx_c) > norminv(1-x),
% with err_bnd = desired_accuracy/n_lin_consts/10.
%
% USAGE: See getUnderapproxStochReachAvoidSet,
% computeCcLowerBoundStochReachAvoidPwlRisk.
%
% =============================================================================
%
% 
% Inputs: None
% -------
%   max_delta        - This is the maximum tolerance for violation of the joint
%                       chance constraint | risk allocation can not exceed this
%                       value
%   desired_accuracy - Accuracy expected from the chance constraint formulation
%                       with n_lin_consts, no. of individual chance constraints
%   n_lin_consts     - No. of individual chance constraints
%   
%
% Outputs:
% --------
%   cdf_approx_m - Secant slopes that will overapproximate norminv(1-x)
%   cdf_approx_c - Secant y-intercepts that will overapproximate norminv(1-x)
%   lb_phiinv    - Lower bound on x in norminv(1-x) for which the provided PWA
%                   approximation
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
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
%
%

    % Max bounds we have tested this for
    compute_phiinv_lb = 1E-8;
    compute_phiinv_ub = 0.44;
    
    % Interested bounds are desired_accuracy/N_{ineq} and max_delta
    lb_phiinv = desired_accuracy/n_lin_consts/10;
    upper_bound_phiinv = max_delta; 

    % Throw an error if the requested bounds are beyond the PWA parameters
    if compute_phiinv_lb > lb_phiinv ||...
        compute_phiinv_ub < upper_bound_phiinv
        exc = SrtInvalidArgsError(['Requested bounds exceed limits.']);
        throwAsCaller(exc);
    end
    
    if ~exist('SReachTools_data.mat','file')
        % Create the look up table for norminv pdf if none found in
        % helperfunctions
        create_PWAapprox_norminv(compute_phiinv_lb, compute_phiinv_ub,...
            desired_accuracy);
    end

    % This matfile will have norminv's cdf_approx_m,cdf_approx_c, and knots
    load('SReachTools_data','norminvcdf_approx_m','norminvcdf_approx_c',...
        'norminv_knots');

    % Given the desired lower bound and upper bound, find what pieces are needed
    % -1 added to knots_ub_indx because each piece is associated with the left
    % hand side knot.
    knots_lb_indx = find(norminv_knots >= lb_phiinv,1);
    knots_ub_indx = find(norminv_knots >= upper_bound_phiinv,1)-1;
    cdf_approx_m = norminvcdf_approx_m(knots_lb_indx:knots_ub_indx);
    cdf_approx_c = norminvcdf_approx_c(knots_lb_indx:knots_ub_indx);
end

function create_PWAapprox_norminv(lb_phiinv,...
    upper_bound_phiinv, maxlierror_phiinv)
    % create_PWAapprox_norminv creates a PWA approximation of normcdfinv

    % phiinv(1-x) definition: https://www.mathworks.com/help/stats/norminv.html
    % And using erfinv instead of erfcinv since we need \Phi^{-1}(1-x)
    phiinv = @(z) sqrt(2)* erfinv(2*(1 - z) -1 );

    % Concave function required: So negate phiinv and it becomes monotone inc
    function_handle = @(z) -phiinv(z);
    fun_hessian_monotone_phiinv = 'mono-inc';

    % get PWA approximation
    [~,~,PWA_negphiinv_underapprox_m, PWA_negphiinv_underapprox_c,...
        norminv_knots] =...
            getPWAOverAndUnderApprox(lb_phiinv,...
                upper_bound_phiinv,...
                maxlierror_phiinv,...
                function_handle,...
                fun_hessian_monotone_phiinv);

    % negate the PWA underapproximation of the concave representation to get the
    % PWA overapproximation of the convex function
    norminvcdf_approx_m = - PWA_negphiinv_underapprox_m';
    norminvcdf_approx_c = - PWA_negphiinv_underapprox_c';

    % Save the data in the helperFunctions folder
    save('../src/helperFunctions/SReachTools_data',...
        'norminvcdf_approx_m','norminvcdf_approx_c','norminv_knots');
end
