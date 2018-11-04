function approx_level_set = SReachSetLag(method_str, sys, prob_thresh, ...
    safety_tube, options)
% Get approximate level set using lagrangian methods
% ============================================================================
%
% This function will get the approximate prob_thresh level set for a stochastic
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
%   method_str  - Lagrangian method,
%                   'lag-over'  -- Lagrangian Overapproximation
%                   'lag-under' -- Lagrangian Underapproximation
%   sys         - LtiSystem object
%   prob_thresh - Probability threshold
%   safety_tube - Tube object
%   options     - Struct of reach set options, see SReachSetOptions
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
    inpar.addRequired('method_str', @(x) any(validatestring(x, ...
        {'lag-under', 'lag-over'})));
    inpar.addRequired('sys', @(x) validateattributes(x, {'LtiSystem', ...
        'LtvSystem'}, {'nonempty'}));
    inpar.addRequired('prob_thresh', @(x) validateattributes(x, {'numeric'}, ...
        {'>=', 0, '<=', 1}));
    inpar.addRequired('safety_tube', @(x) validateattributes(x, {'Tube'}, ...
        {'nonempty'}));

    try
        inpar.parse(method_str, sys, prob_thresh, safety_tube);
    catch err
        exc = SrtInvalidArgsError.withFunctionName();
        exc = exc.addCause(err);
        throwAsCaller(exc);
    end
    
    % 1. Ensure sys is Gaussian-perturbed system
    % 2. Options string match
    otherInputHandling(method_str, sys, options);
    
    % Obtain time horizon
    time_horizon = length(safety_tube) - 1;
    
    if time_horizon > 0 && prob_thresh > 0        
        switch(lower(method_str))
            case 'lag-under'
                % get bounded scaled_disturbance set
                bounded_set = SReachSetLagBset(sys, ...
                    prob_thresh^(1/time_horizon), options);
                
                % get underapproximated level set (robust effective target)
                approx_level_set = getSReachLagUnderapprox(sys, ...
                    safety_tube, bounded_set);
            case 'lag-over'
                % get bounded disturbance set
                bounded_set = SReachSetLagBset(sys, ...
                    (1-prob_thresh)^(1/time_horizon), options);
                % get overapproximated level set (augmented effective target)
                approx_level_set = getSReachLagOverapprox(sys, ...
                    safety_tube, bounded_set);
            otherwise
                throw(SrtInvalidArgsError('Unhandled method_str: %s', ...
                    method_str));
        end
    elseif time_horizon == 0 || prob_thresh == 0
        % return set in target tube if length is 1
        approx_level_set = safety_tube(1);
    else
        throw(SrtDevError('Unknown problem configuration'));
    end
end

function otherInputHandling(method_str, sys, options)
    
    % Ensure Gaussian-perturbed system
    validateattributes(sys.dist, {'RandomVector'}, {'nonempty'});
    validatestring(sys.dist.type, {'Gaussian'}, {'nonempty'});
    
    % Check if prob_str and method_str are consistent        
    if ~strcmpi(options.prob_str,'term')
        throwAsCaller(...
            SrtInvalidArgsError('Mismatch in prob_str in the options'));
    end
    if ~strcmpi(options.method_str,method_str)
        throwAsCaller(...
            SrtInvalidArgsError('Mismatch in method_str in the options'));
    end        
end
    
