function varargout = SReachSetLag(method_str, sys, prob_thresh, safety_tube,...
    options)
% Get approximate stochastic reach set using Lagrangian methods
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
%   approx_set - Polyhedron object for the over-/under-approximation of the 
%                stochastic reach set
%   approx_tube- [Optional] Tube comprising of an over-/under-approximation of
%                the stochastic reach sets across the time horizon
%   bounded_set- [Optional] Bounded disturbance set which was used to
%                robustify the computation (lag-under) or 
%                augment the input set (lag-over)
%
% Notes:
% ------
% * compute_style of `support` and method of `lag-over` will return only the 
%   approx_set.
% * While 'Gaussian' disturbance can have options.bound_set_method be 'polytope'
%   or 'ellipsoid', 'UserDefined' disturbance requires options.bound_set_method
%   to be 'polytope'.
% ============================================================================
% 
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://sreachtools.github.io/license/
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
    
    if time_horizon == 0 || prob_thresh == 0
        % return set in target tube if length is 1
        varargout{1} = safety_tube(1);
        varargout{2} = safety_tube;
        varargout{3} = [];
    elseif time_horizon > 0 && prob_thresh > 0   
        if options.verbose
            under_over=strsplit(method_str,'-');
            fprintf('Computing Lagragian %s approximation\n\n', under_over{2});
        end
        switch(lower(method_str))
            case 'lag-under'
                % get bounded scaled_disturbance set
                bounded_set = SReachSetLagBset(sys, ...
                    prob_thresh^(1/time_horizon), options);
                
                [approx_set, approx_tube] = getSReachLagUnderapprox(...
                    sys, safety_tube, bounded_set, options);
                varargout{1} = approx_set;
                varargout{2} = approx_tube;
                varargout{3} = bounded_set;
            case 'lag-over'
                % get bounded disturbance set
                bounded_set = SReachSetLagBset(sys, ...
                    (1-prob_thresh)^(1/time_horizon), options);
                % get overapproximated set (augmented effective target)
                switch lower(options.compute_style)
                    case 'vfmethod'
                        [approx_set, approx_tube] = getSReachLagOverapprox(...
                            sys, safety_tube, bounded_set, options);
                        varargout{1} = approx_set;
                        varargout{2} = approx_tube;
                        varargout{3} = bounded_set;
                    case 'support'
                        if nargout >= 2
                            throw(SrtInvalidArgsError(['Too many output ',...
                                'arguments.\ncompute_style = support can ',...
                                'not compute overapproximation REACH TUBE.']));
                        end
                        approx_set = getSReachLagOverapprox(sys, safety_tube,...
                            bounded_set, options);
                        varargout{1} = approx_set;
                        varargout{2} = [];
                        varargout{3} = bounded_set;
                end                
            otherwise
                throw(SrtInvalidArgsError('Unhandled method_str: %s', ...
                    method_str));
        end
    else
        throw(SrtDevError('Unknown problem configuration'));
    end
end

function otherInputHandling(method_str, sys, options)
    
    % Ensure stochastic system
    validateattributes(sys.dist, {'RandomVector'}, {'nonempty'},...
        'SReachSetLag/otherInputHandling', 'sys.dist');
    
    % Check if prob_str and method_str are consistent        
    if ~strcmpi(options.prob_str,'term')
        throwAsCaller(...
            SrtInvalidArgsError('Mismatch in prob_str in the options'));
    end
    if ~strcmpi(options.method_str,method_str)
        throwAsCaller(...
            SrtInvalidArgsError('Mismatch in method_str in the options'));
    end        

    % Based on the RandomVector, decide allowable methods
    validatestring(sys.dist.type, {'Gaussian','UserDefined'}, {'nonempty'},...
        'SReachSetLag/otherInputHandling', 'sys.dist.type');
end
    
