function [stoch_reach_prob, opt_input_vec, opt_input_gain, varargout] =...
    SReachPoint(prob_str, method_str, sys, initial_state, safety_tube, varargin)
% Solve the stochastic (first/terminal) reach problem approximately from a given
% initial state using a host of techniques
% =============================================================================
%
% SReachPoint computes an approximation to the first/terminal stochastic reach
% problems. 
%
% TODO: Talk about first and terminal hitting time problems
%
% This function is a compilation of various techniques proposed in the
% literature:
%
% 1. Convex chance-constrained-based approach (chance-open):
%
%    High-level desc.   : Use Boole's inequality, Gaussian random vector, and
%                         piecewise linear approximation of the inverse of the
%                         standard normal cumulative density function to create
%                         a linear program-based approximation to the original
%                         optimization
%    Approximation      : Guaranteed underapproximation
%    Controller type    : Open-loop controller that satisfies the hard
%                         input bounds 
%    Optimality         : Optimal open-loop controller for the
%                         underapproximation problem due to convexity guarantees
%    SReachTool function: getCcOpenPoint
%    Dependency (EXT)   : CVX
%    Dependency (MATLAB): Symbolic toolbox
%    Paper              : a. Lesser, Oishi, Erwin TODO.
%                         b. A. Vinod and M. Oishi, HSCC 2018 TODO
%
% 2. Convex chance-constrained-based approach (chance-affine):
%
%    High-level desc.   : Use Boole's inequality, Gaussian random vector,
%                         hyperbolic constraints-to-second order cone constraint
%                         reformulation, and piecewise linear approximation of
%                         the inverse of the standard normal cumulative density
%                         function to create a second-order cone program-based
%                         difference-of-convex optimization problem
%    Controller type    : A history-dependent affine controller that satisfies
%                         softened input constraints (controller satisfies the
%                         hard input bounds upto a user-specified probabilistic
%                         threshold)
%    Optimality         : Suboptimal affine controller for the
%                         underapproximation problem due to the use of
%                         difference-of-convex
%    Approximation      : Guaranteed underapproximation
%    SReachTool function: getCcClosedPoint
%    Dependency (EXT)   : CVX
%    Dependency (MATLAB): Symbolic toolbox
%    Paper              : A. Vinod and M. Oishi, HSCC 2018 TODO
%
% 3. Fourier transform + Patternsearch (genzps-open):
%
%    High-level desc.   : Maximize the multivariate Gaussian integral over a
%                         polytope, evaluated using Genz's algorithm, and
%                         optimize the nonlinear (log-concave) problem using
%                         MATLAB's patternsearch
%    Approximation      : Approximate upto a user-specified tolerance
%    Controller type    : Open-loop controller that satisfies the hard input
%                         bounds
%    Optimality         : Optimal open-loop controller for the
%                         underapproximation problem due to convexity guarantees
%    Dependency (MATLAB): Global Optimization toolbox (for patternsearch)
%    SReachTool function: getFtOpenPoint
%    Paper              : A. Vinod and M. Oishi, "Scalable Underapproximation
%                         for Stochastic Reach-Avoid Problem for
%                         High-Dimensional LTI Systems using Fourier
%                         Transforms," in IEEE Control Systems Letters, 2017.
%
% 4. Scenario-based approach (scenario-open):
%
%    High-level desc.   : Sample scenarios based on the additive noise and solve
%                         a mixed-integer linear program to make the maximum
%                         number of scenarios satisfy the reachability objective
%    Approximation      : No direct approximation guarantees. Accuracy improves
%                         as the number of scenarios considered increases.
%    Controller type    : Open-loop controller that satisfies the hard input
%                         bounds
%    Optimality         : Optimal (w.r.t scenarios drawn) open-loop controller
%                         for the underapproximation problem 
%    Dependency (EXT)   : CVX
%    SReachTool function: getScenarioPoint
%    Paper              : Lesser, Oishi, Erwin TODO.
%
%
% USAGE: TODO
%
% =============================================================================
%
% [stoch_reach_prob, opt_controller, varargout] = SReachPoint(prob_str,...
%    method_str, sys, init_state, safety_tube, options)
% 
% Inputs:
% -------
%   prob_str    - String specifying the problem of interest. For each case, we
%                 compute the optimal value function that maps initial states
%                 to different maximal reach probabilities
%                     1. 'first' : Stay within the safety_tube and reach the
%                                  target set early if possible
%                     2. 'term' : Stay within the safety_tube
%   method_str  - Solution technique to be used.
%                     'chance-open'  -- Convex chance-constrained approach for
%                                       an open-loop controller synthesis
%                     'chance-affine'-- Convex chance-constrained approach for
%                                       an affine controller synthesis
%                     'genzps-open'  -- Genz's algorithm + Patternsearch
%                     'scenario-open'-- Scenario-based 
%   sys         - System description (LtvSystem/LtiSystem object)
%   init_state  - Initial state for which the maximal reach probability must be
%                 evaluated (A numeric vector of dimension sys.state_dim)
%   safety_tube - Collection of (potentially time-varying) safe sets that define
%                 the safe states (TargetTube object)
%   options     - Collection of user-specified options for each of the solution
%                 (Matlab struct created using SReachOptions)
%
% Outputs:
% --------
%   stoch_reach_prob 
%               - Approximation (underapproximation, in some cases) to the
%                 first/terminal stochastic reach problem
%   opt_controller 
%               - Controller U=MW+d for a concatenated input vector 
%                   U = [u_0; u_1; ...; u_{N-1}] and concatenated disturbance
%                   vector W=[w_0; w_1; ...; w_{N-1}]. 
%                 MATLAB struct with the following members:
%                   - M: Affine controller gain matrix of dimension
%                       (sys.input_dim*N) x (sys.dist_dim*N)
%                   - d: Open-loop controller: column vector dimension
%                       (sys.input_dim*N) x 1
%                 The feedback gain matrix M is set to [] for methods that look
%                 for open-loop controllers. 
%   TODO: Optional outputs
%
% Notes:
% * See @LtiSystem/getConcatMats for more information about the
%   notation used.
% * If an open_loop policy is desired arranged in increasing time columnwise,
%   use the following command
%   TODO
%   >>> optimal_open_loop_control_policy = reshape(opt_controller,...
%           sys.input_dim, time_horizon);
% 
% =============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
%
%

    % Input parsing
    valid_prob = {'first','term'};
    valid_method= {'chance-open','chance-affine','genzps-open','scenario-open'};

    inpar = inputParser();
    inpar.addRequired('prob_str', @(x) any(validatestring(x,valid_prob)));
    inpar.addRequired('method_str', @(x) any(validatestring(x,valid_method)));
    inpar.addRequired('sys', @(x) validateattributes(x,...
        {'LtiSystem','LtvSystem'}, {'nonempty'}));
    inpar.addRequired('init_state', @(x) validateattributes(x, {'numeric'},...
        {'vector'}));
    inpar.addRequired('safety_tube',@(x) validateattributes(x,{'TargetTube'},...
        {'nonempty'}));

    try
        inpar.parse(prob_str, method_str, sys, initial_state, safety_tube);
    catch err
        exc = SrtInvalidArgsError.withFunctionName();
        exc = exc.addCause(err);
        throwAsCaller(exc);
    end
        
    % Check if safe set contains the initial state
    if ~safety_tube(1).contains(initial_state)
        % Stochastic reach-avoid probability is zero and no admissible open-loop
        % policy exists, if given an unsafe initial state
        stoch_reach_prob = 0;
        opt_input_vec = nan(sys.input_dim * time_horizon, 1);
        opt_input_gain = [];        
    else
        % Ensure that options are provided are appropriate
        options = otherInputHandling(prob_str,method_str, varargin);
        
        % Dependig on method_str call the appropriate solution technique
        switch(lower(method_str))
            case 'genzps-open'
                % Patternsearch and Fourier transform-based open-loop
                % underapproximation of the stochastic reach-avoid problem
                [stoch_reach_prob, opt_input_vec] = SReachPointGpO(sys,...
                    initial_state, safety_tube, options.desired_accuracy,...
                    options.PSoptions);
                 opt_input_gain = [];
            case 'chance-open'
                % Chance-constrained formulation with piecewise-linear 
                % approximations to compute open-loop controller (LP)
                [stoch_reach_prob, opt_input_vec, risk_alloc_state] =...
                    SReachPointCcO(sys, initial_state, safety_tube,...
                        options.pwa_accuracy);
                 opt_input_gain = [];
                 varargout{1} = risk_alloc_state;
            case 'scenario-open'
                % Chance-constrained formulation with piecewise-linear 
                % approximations to compute open-loop controller (LP)
                [stoch_reach_prob, opt_input_vec] = SReachPointCcO(sys,...
                    initial_state, safety_tube, options.desired_accuracy);
                 opt_input_gain = [];
            case 'chance-affine'
                % Chance-constrained formulation with piecewise-linear 
                % approximations to compute affine-loop controller (SOC program)
                [stoch_reach_prob_affine, opt_input_vec, opt_input_gain,...
                    risk_alloc_state, risk_alloc_input] = SReachPointCcA(sys,...
                        initial_state, safety_tube, options);
                stoch_reach_prob = stoch_reach_prob_affine; %TODO
                varargout{1} = risk_alloc_state;
                varargout{2} = risk_alloc_input;
        end
        if stoch_reach_prob<0
            stoch_reach_prob = 0;
            opt_input_vec = nan(sys.input_dim * time_horizon, 1);
            opt_input_gain = [];
            if strcmpi(method_str,'chance-open') ||...
                    strcmpi(method_str,'chance-affine')    %TODO
                warn_str = ['''chance-open/chance-affine'' works only',...
                    ' if the maximal reach probability is above 0.5.'];
            else
                warn_str = 'The problem may be infeasible';
            end
            % Need to add sprintf (despite MATLAB's editor warning) for new line
            warning(sprintf(method_str,' returned a trivial lower bound.\n',...
                warn_str));             
        end
    end
end

function options = otherInputHandling(prob_str,method_str, additional_args)
    % input handling for options

    % Check if options provided is appropriate or generate
    if isempty(additional_args)
        % No options provided! Create one.
        options = SReachPointOptions(prob_str,method_str);
    elseif length(additional_args)==1
        % Options provided! Check if prob_str and method_str are consistent
        options = additional_args{1};
        if ~strcmpi(options.prob_str,prob_str)
            throwAsCaller(SrtInvalidArgsError('Mismatch in prob_str in the options'));
        end
        if ~strcmpi(options.method_str,method_str)
            throwAsCaller(SrtInvalidArgsError('Mismatch in method_str in the options'));
        end        
    else
        % Lots more than needed!
        throwAsCaller(SrtInvalidArgsError('Too many input arguments'));
    end
end