function approx_level_set = getApproxStochasticLevelSetViaLagrangian(sys, ...
    beta, target_tube, approx_type, method, varargin)
% SReachTools/stochasticReachAvoid/getApproxStochasticLevelSetViaLagrangian: 
% Get approximate level set using lagrangian methods
% ============================================================================
%
% This function will get the approximate beta level set for a stochastic
% discrete time system using the Lagrangian methods in 
%     J. D. Gleason, A. P. Vinod, M. M. K. Oishi, "Underapproximation of 
%     Reach-Avoid Sets for Discrete-Time Stochastic Systems via Lagrangian 
%     Methods," in Proceedings of the IEEE Conference on Decision and Control, 
%     2017
%
% Usage: see examples/doubleIntegratorLevelSetApprox.m
%        or  examples/lagrangianApproximations.m
% 
% ============================================================================
%
% Inputs:
% -------
%   sys              - LtiSystem object
%   beta             - Probability threshold
%   target_tube      - Cell array of polyhedron objects
%   approx_type      - Approximation type, either 'overapproximation' or
%                      'underapproximation'
%   method, varargin - See help getBoundedSetForDisturbance
%
% Outputs:
% --------
%   approx_level_set - Polyhedron object
%
% ============================================================================
% 
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 
% 

    % verify inputs
    inpar = inputParser();
    inpar.addRequired('sys', @(x) validateattributes(x, {'LtiSystem',...
        'LtvSystem'}, {'nonempty'}));
    inpar.addRequired('beta', @(x) validateattributes(x, {'numeric'}, ...
        {'>=', 0, '<=', 1}));
    inpar.addRequired('target_tube', ...
        @(x) validateattributes(x, {'TargetTube'}, {'nonempty'}));
    inpar.addRequired('approx_type', @(x) any(validatestring(x, ...
        {'underapproximation', 'overapproximation'})));
    inpar.addRequired('method', @(x) any(validatestring(x, ...
        {'random', 'box', 'load'})));

    inpar.parse(sys, beta, target_tube, approx_type, method);

    % additional non-input parser validations
    validateattributes(sys.dist, {'RandomVector','StochasticDisturbance'}, ...
        {'nonempty'});
    
    % validate that all elements of the target_tube are polyhedron
    for i = 1:length(target_tube)
        validateattributes(target_tube(i), {'Polyhedron'}, {'nonempty'});
    end
    
    validateattributes(approx_type, {'char', 'string'}, {'nonempty'});
    switch(approx_type)
        case 'underapproximation'
            do_underapprox = true;
        case 'overapproximation'
            do_underapprox = false;
        otherwise
            error('SReachTools:invalidArgs', ['Input ''approx_type'' must be ', ...
                'either ''underapproximation'' or ''overapproximation'', ', ...
                'see help getApproxStochasticLevelSetViaLagrangian']);
    end
    
    validateattributes(method, {'char', 'string'}, {'nonempty'});
    
    if length(target_tube) > 1
        if do_underapprox
            % Perform underapproximation

            % get bounded disturbance set
            bounded_set = getBoundedSetForDisturbance(sys.dist, ...
                length(target_tube)-1, beta, method, varargin{:});

            % get underapproximated level set (robust effective target)
            approx_level_set = getRobustEffTarget(sys, target_tube, ...
                bounded_set);
        else
            % get bounded disturbance set
            bounded_set = getBoundedSetForDisturbance(sys.dist, ...
                length(target_tube)-1, (1-beta), method, varargin{:});
            
            % get overapproximated level set (augmented effective target)
            approx_level_set = getAugEffTarget(sys, target_tube, ...
                bounded_set);
        end
    else
        % return set in target tube if length is 1
        approx_level_set = target_tube{1};
    end
end
