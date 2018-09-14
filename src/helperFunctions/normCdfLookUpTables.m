function normCdfLookUpTables(refresh_tables)
% Computes lookup tables for Gaussian quantile functions
% ============================================================================
% 
% Computes piecewise affine overapproximation (lookup tables) for the inverse of
% the standard normal cumulative density function (quantile function).
% It uses getPWAOverAndUnderApprox to create this overapproximation upto various
% error margins. See the following paper for more details.
%
% Abraham P. Vinod, Vignesh Sivaramakrishnan, and Meeko M. K. Oishi,
% "Piecewise-Affine Approximation-Based Risk Allocation for Gaussian Joint
% Chance Constraints", in American Control Conference, 2019 (submitted). TODO
%
% ============================================================================
% 
% normCdfLookUpTables(refresh_tables)
%
% Inputs:
% -------
% refresh_tables  - Set to 1 to clear existing mat files and create new tables
%
% Outputs:
% --------
%
% Notes:
% ------
% * MATLAB DEPENDENCY: Uses MATLAB's Symbolic toolbox for ease in implementation
% * Stores the lookup tables in src/helperFunction
%
% ============================================================================
%
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
%
%
    
    lb_phiinv = 1e-8;
    ub_phiinv = 0.5;
    err_ranges = [1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2];
    nocreate = 1;
    
    if refresh_tables
        delete([getSReachToolsHome() 'src/helperFunctions/SReachTools_norminvcdf_*.mat']);
    end
    SReachToolsHome = getSReachToolsHome();
    for maxerror = err_ranges
        mat_identifier = replace(replace(num2str(maxerror,'%1.0e'),'.','x'),'-','_');
        mat_str = [SReachToolsHome 'src/helperFunctions/SReachTools_norminvcdf_' mat_identifier '.mat'];
    
        if ~exist(mat_str,'file')    
            if nocreate
                disp(['Creating piecewise affine overapproximations for norminv(1-x) for x in [',...
                num2str(lb_phiinv) ',' num2str(ub_phiinv) ']. This is a one-time computation.']);
                nocreate = 0;
            end
            create_PWAapprox_norminv(lb_phiinv, ub_phiinv, maxerror);        
        end
    end
    if nocreate
        disp('Found lookup tables for norminv! Skipped creating new tables.');
    end
end


function create_PWAapprox_norminv(lb_phiinv, ub_phiinv, maxlierror_phiinv)
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
                ub_phiinv,...
                maxlierror_phiinv,...
                function_handle,...
                fun_hessian_monotone_phiinv);

    % negate the PWA underapproximation of the concave representation to get the
    % PWA overapproximation of the convex function
    norminvcdf_approx_m = - PWA_negphiinv_underapprox_m';
    norminvcdf_approx_c = - PWA_negphiinv_underapprox_c';
    fprintf('For max error: %1.2e, we needed %d affine pieces\n', maxlierror_phiinv,length(norminvcdf_approx_m));
    % Save the data in the helperFunctions folder
    mat_identifier = replace(replace(num2str(maxlierror_phiinv,'%1.0e'),'.','x'),'-','_');
    save([getSReachToolsHome() 'src/helperFunctions/SReachTools_norminvcdf_' mat_identifier '.mat'],...
        'norminvcdf_approx_m','norminvcdf_approx_c','norminv_knots');
end
