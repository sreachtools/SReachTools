function [cdf_approx_m, cdf_approx_c, lb_phiinv, useful_knots] =...
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
%   useful_knots - Breakpoints of the PWA overapproximation, i.e., points
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
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
%
%
    
    % Max bounds we have tested this for
    compute_phiinv_lb = 1E-8;
    compute_phiinv_ub = 0.5;
    err_ranges = [1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2];

    % Pick the largest lower bound to the desired accuracy
    err_indx = find(desired_accuracy>=err_ranges,1,'last');
    if isempty(err_indx)
        warning('Requested accuracy of %1.2e, but using %1.0e',desired_accuracy, err_ranges(1));
        err_indx = 1;
    end
    
    % Interested bounds are desired_accuracy/N_{ineq} and max_delta
    lb_phiinv = desired_accuracy/n_lin_consts/10;
    upper_bound_phiinv = max_delta;
    
    % Throw an error if the requested bounds are beyond the PWA parameters
    if compute_phiinv_lb > lb_phiinv || compute_phiinv_ub < upper_bound_phiinv
        exc = SrtInvalidArgsError(['Requested bounds exceed limits.']);
        throw(exc);
    end
    
    mat_identifier = replace(replace(num2str(err_ranges(err_indx),'%1.0e'),'.','x'),'-','_');
    mat_str = [getSReachToolsHome() 'src/helperFunctions/SReachTools_norminvcdf_' mat_identifier '.mat'];
    if exist(mat_str,'file')
        % This matfile will have norminv's cdf_approx_m,cdf_approx_c, and knots
        load(mat_str,'norminvcdf_approx_m','norminvcdf_approx_c', 'norminv_knots');
    else
        normCdfLookUpTables(0);
        % This matfile will have norminv's cdf_approx_m,cdf_approx_c, and knots
        load(mat_str,'norminvcdf_approx_m','norminvcdf_approx_c', 'norminv_knots');
    end
    % Given the desired lower bound and upper bound, find what pieces are needed
    % -1 added to knots_ub_indx because each piece is associated with the left
    % hand side knot.
    knots_lb_indx = find(norminv_knots >= lb_phiinv,1);
    knots_ub_indx = find(norminv_knots >= upper_bound_phiinv,1)-1;
    cdf_approx_m = norminvcdf_approx_m(knots_lb_indx:knots_ub_indx);
    cdf_approx_c = norminvcdf_approx_c(knots_lb_indx:knots_ub_indx);
    useful_knots = norminv_knots(knots_lb_indx:knots_ub_indx);
end


