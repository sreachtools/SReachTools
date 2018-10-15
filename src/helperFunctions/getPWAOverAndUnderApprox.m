function [PWA_overapprox_m,...
          PWA_overapprox_c,...
          PWA_underapprox_m,...
          PWA_underapprox_c,...
          knots_underapprox] = getPWAOverAndUnderApprox(lb,...
            ub,...
            desired_accuracy,...
            g_matlabfun,...
            hessian_monotone)
% Computes a PWA over- and underapproximation of a concave function
% ============================================================================
% 
% This function implements the piecewise affine (PWA) approximation algorithm
% presented in
%
% Abraham P. Vinod, Vignesh Sivaramakrishnan, and Meeko M. K. Oishi,
% "Piecewise-Affine Approximation-Based Risk Allocation for Gaussian Joint
% Chance Constraints", in American Control Conference, 2019 (submitted). TODO
%
% It uses linear Lagrange interpolation formula for the underapproximation and
% the first-order Taylor series for the overapproximation computation. 
%
% Usage: See computeNormCdfInvOverApprox.m
%
% ============================================================================
% 
% [PWA_overapprox_m,...
%  PWA_overapprox_c,...
%  PWA_underapprox_m,...
%  PWA_underapprox_c,...
%  knots_underapprox] = getPWAOverAndUnderApprox(lb,...
%    ub,...
%    desired_accuracy,...
%    g_matlabfun,...
%    hessian_monotone)
%
% Inputs:
% -------
% lb                - Lower bound on the domain
% ub                - Upper bound on the domain
% desired_accuracy  - Accuracy of the PWA approximation
% g_matlabfun       - Function handle for the concave function that is to be
%                       approximated
% hessian_monotone  - Montonicity of Hessian ('mono-inc'/'mono-dec')
%
% Outputs:
% --------
% PWA_overapprox_m, PWA_overapprox_c
%                   - Slope and intercept of the PWA overapproximation
% PWA_underapprox_m, PWA_underapprox_c
%                   - Slope and intercept of the PWA underapproximation
% knots_underapprox - Intercept points for PWA underapproximation and given
%                       function
%
% Notes:
% ------
% * MATLAB DEPENDENCY: Uses MATLAB's Symbolic toolbox for ease in implementation
% * Requires the concave function to have a monotone hessian. This enables easy
%       implementation of the error bound via linear Lagrange interpolation
% * Can be used for convex functions as well by negating the input function and
%   negating the output affine functions
% * The overapproximation part of the code assumes a gradient at the upper bound
%   if the curvature is not sufficient enough for the mean value theorem to
%   guarantee a point whose gradient matches the underapproximation slope.
% * This function is implicitly covered in testing in
%   computeNormCdfInvOverApproxTest.m
%
% ============================================================================
%
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
%
%
    
    %% Initialization
    knots_underapprox(1) = lb;
    j = 1;
    fzero_options = optimset('Display','off'); % show iterations
            
    % Ensure that symbolic math toolbox is installed
    v = ver;
    has_syms = any(strcmp(cellstr(char(v.Name)), 'Symbolic Math Toolbox'));
    if ~has_syms
        exc = SrtSetupError(['This function needs MATLAB''s ', ...
            'Symbolic Math Toolbox.']);
        throw(exc);
    end

    %% Symbolic function definitions
    syms x;
    syms h;
    syms hmax;

    % Function definition for the input function
    g = g_matlabfun(x);
    % First derivative
    g1diff=diff(g,1);
    % Second derivative
    g2diff=diff(g,2);
    % Due to concavity of g, the error becomes more negative for larger h
    
    % Iterate till we reach the end
    while knots_underapprox(j)< ub
        % Solve h^2 min_{y\in[x_j,x_j+h]} d^2/dh^2 g(y) == -8 * desired_accuracy
        if strcmpi(hessian_monotone,'mono-inc')
            % If monotone increasing,
            % min_{y\in[x_j,x_j+h]} d^2/dh^2 g(y) = d^2/dh^2 g(x_j)
            % Solve for h very easily 
            % Added abs to allow for x^2
            hval = sqrt(abs(-8*desired_accuracy/...
                double(subs(g2diff,x,knots_underapprox(j)))));
            knots_underapprox(j+1) = min([knots_underapprox(j) + hval,ub]);
        elseif strcmpi(hessian_monotone,'mono-dec')
            % If monotone decreasing, 
            % min_{y\in[x_j,x_j+h]} d^2/dh^2 g(y) = d^2/dh^2 g(x_j + h)
            % So, slightly harder root-finding problem
            g_maxerror=subs(g2diff,x+h)*h^2 + 8*desired_accuracy;
            g_maxerror_at_j = subs(g_maxerror,x,knots_underapprox(j));
            g_maxerror_at_j_mf = @(z) double(subs(g_maxerror_at_j,z));
            % Set up the search interval to not exceed ubdelta
            search_interval = [0 ub-knots_underapprox(j)];
            % Solve for h
            try
                hval = fzero(g_maxerror_at_j_mf,search_interval,fzero_options);
                knots_underapprox(j+1) = knots_underapprox(j) + hval;
            catch
                % If errored, then the search interval doesn't see a sign flip
                % => we have reached the end
                knots_underapprox(j+1) = ub;
            end
        else
            exc = SrtInvalidArgsError(['Concave function needs to have ',...
                'monotone (increasing/decreasing) hessian']);
            throwAsCaller(exc);
        end
        % disp(knots_underapprox(j));
        
        %% Construction of the underapproximation via Lagrange interpolation
        y_2 = double(subs(g,knots_underapprox(j+1)));
        y_1 = double(subs(g,knots_underapprox(j)));
        x_2 = knots_underapprox(j+1);
        x_1 = knots_underapprox(j);

        % Lagrange linear interpolation
        PWA_underapprox_m(j) = (y_2 - y_1)/(x_2 - x_1);
        PWA_underapprox_c(j) = y_1 - PWA_underapprox_m(j) * x_1;
        
        %% Construction of the overapproximation via first-order Taylor series
        % Set up the gradient function
        g1diff_at_j = g1diff - PWA_underapprox_m(j);
        % Using matlabFunction directly is not advised 
        % TODO: Write blogpost
        g1diff_at_j_mf = @(z) double(subs(g1diff_at_j,z));
        % Set up the search interval as [x(j), x(j+1)]
        search_interval = [x_1  x_2];
        % Search for hvalgrad such that f'(x(j)+hvalgrad) = c_{j} ---
        % existence guaranteed by mean value theorem
        try
            [x_grad_match] = fzero(g1diff_at_j_mf,search_interval,fzero_options);
            %cdf_underapprox_m(j) is the same as double(subs(g1diff,x_grad_match));
            PWA_overapprox_m(j) = PWA_underapprox_m(j);
            PWA_overapprox_c(j) = double(subs(g,x_grad_match)) -...
                PWA_overapprox_m(j) * x_grad_match;        
        catch
            if x_2 < ub
                exc = SrtInternalError('Was expecting the endpoint! Internal error?!');
                throw(exc)
            else
                % Towards the endpoint, the curve may not accommodate MVT
                PWA_overapprox_m(j) = double(subs(g1diff,x_2));
                PWA_overapprox_c(j) = double(subs(g,x_2)) -...
                    PWA_overapprox_m(j) * x_2;        
            end
        end        
        %% Increment j
        j = j+1;       
    end
end
