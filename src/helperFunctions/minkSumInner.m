function polyInner = minkSumInner(polyA, polyB, varargin)    
    % Computes the scaling and shifting for the template polytope polyC
    % such that it is completely contained in the polytope polyA + polyB
    % If no template is provided, we use polyA as the template

    if nargin > 3
        throwAsCaller(SrtInvalidArgsError('Too many input arguments'));
    elseif nargin == 3
        polyC = varargin{1};
    else
        polyC = polyA;
    end
    polyA.minHRep();
    polyB.minHRep();
    polyC.minHRep();
    if ~isempty(polyC.He)
        throwAsCaller(SrtInvalidArgsError(['Template must have only ',...
            'inequalities']));
    end
    if isEmptySet(polyA) || isEmptySet(polyB) || isEmptySet(polyC)
        throwAsCaller(SrtInvalidArgsError('Empty set provided as arguments'));
    end
    if ~((polyA.Dim == polyB.Dim) && (polyB.Dim == polyC.Dim))
        throwAsCaller(SrtInvalidArgsError('Mismatch in polytope dimensions'));
    end
    
    % To avoid the wierd CVX error of undefined 'variable'
    cvx_clear
    
    % Ensure that polyC is "centered" (This way uniform scaling is not partial)
    polyC = polyC - polyC.chebyCenter.x;    
    
    % To simplify the problem, we will not optimize for the shift but
    % simply use the minkowski sum of Chebyshev centers (deepest points)
    % TODO: Can we guarantee that Minkowski sum's chebyshev center is this
    % point exactly?
    xc_guess = polyA.chebyCenter().x + polyB.chebyCenter().x;
    
    
    % Bisect
    t_ub = 10;
    t_lb = 0;
    bisect_tol = 5e-2;
    myzero_lb = -1e-6;
    
    contain_chk = @(t) outer_opt(t, xc_guess, polyA, polyB, polyC) >= myzero_lb;
%     contain_chk = @(t) bilin_opt(t, xc_guess, polyA, polyB, polyC) >= myzero_lb;
    
    if contain_chk(t_ub)
        warning('SReachTools:runtime',['10x magnification of the template',...
            ' is a subset of the Minkowski sum. Please use a larger ',...
            'polytope for better results.']);
        t_opt = t_ub;
    else
        fprintf('Bisection upper bound of t=%1.3f is valid\n',t_ub);
        t_opt = 0;
        while abs(t_ub - t_lb) > bisect_tol
            t_try = (t_ub + t_lb)/2;
            if contain_chk(t_try)
                t_opt = t_try;
                t_lb = t_try;
                fprintf('t=%1.3f:   Feasible\n',t_try);
            else
                t_ub = t_try;
                fprintf('t=%1.3f: Infeasible\n',t_try);
            end
        end
    end    
    if t_opt < bisect_tol
        warning('SReachTools:runtime', sprintf(['%1.2fx magnification of ',...
            'the template is also not a subset of the Minkowski sum. Please',...
            ' use a smaller polytope for a non-trivial results.'],bisect_tol));
    end
    xc_opt = xc_guess;
    polyInner = t_opt * polyC + xc_opt;
end

function optval = bilin_opt(t, xc, polyA, polyB, polyC)
    ell = sdpvar(polyA.Dim,1);
    lambda1 = sdpvar(size(polyA.A,1),1);
    lambda2 = sdpvar(size(polyB.A,1),1);
    x = sdpvar(polyC.Dim,1);
    obj_val = sdpvar;
    constraints = [polyA.A'*lambda1 == ell, polyB.A'* lambda2 == ell,...
        lambda1 >= 0, lambda2 >= 0, polyC.A * x <= polyC.b,...
        polyA.b'*lambda1 + polyB.b'*lambda2 - ell'*x*t - ell'*xc <= obj_val,...
        ell<= 1, ell >= -1];
    options = sdpsettings('verbose',1,'solver','bmibnb');
    sol = optimize(constraints, obj_val, options);
    if sol.problem == 3
        % Infeasible due to termination limit
        optval = Inf;
    elseif sol.problem == 2
        % Unbounded, so set to -Inf
        optval = -Inf;
    else
        optval = value(obj_val);
    end
end

function optval = outer_opt(t, xc, polyA, polyB, polyC)
    % Compute the support function vector that can violate containment,
    % i.e. inf_l sup_{lambda, y, z} containment_obj
    % We need this value to be positive for guranteeing containment
    options = optimoptions('fmincon','Display','off');
    [~,optval] = fmincon(@(ell) inner_opt(ell, t, xc, polyA, polyB, polyC),...
        0.5*ones(polyA.Dim,1), [], [], [], [],...
        -ones(polyA.Dim,1), ones(polyA.Dim,1), [], options);
end

function optval = inner_opt(ell, t, xc, polyA, polyB, polyC)
    % Solve the inner optimization problem that ensures containment for a
    % given support function vector ell, scaling t, shift xc and polytopes
    % polyA, polyB (Minkowski summands) and polyC (template polytope that
    % is to be scaled and shifted)    
    cvx_clear
    cvx_begin quiet
        cvx_precision best
        variable lambda(size(polyC.H,1),1);
        variable y(polyA.Dim,1);
        variable z(polyB.Dim,1);
        
        maximize (ell' * (y+z-xc) - polyC.b'*lambda)
        subject to
            polyA.A * y <= polyA.b;
            polyA.Ae * y == polyA.be;
            polyB.A * z <= polyB.b;
            polyB.Ae * z == polyB.be;
            lambda >= 0;
            polyC.A' * lambda == t * ell;            
    cvx_end    
    optval = cvx_optval;  
    if strcmpi(cvx_status,'Solved')        
    elseif strcmpi(cvx_status,'Unbounded')
        % Setting it to a large pos value so that fmincon won't freak out
        optval = 1e10;
    elseif strcmpi(cvx_status,'Infeasible')
        % Setting it to a large neg value so that fmincon won't freak out
        optval = -1e10;
    else
        keyboard
    end
end