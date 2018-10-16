function [approx_reach_prob, opt_input_vec, opt_input_gain, varargout] =...
    SReachPoint(prob_str, method_str, sys, initial_state, safety_tube, varargin)
% Solve the stochastic (first/terminal) reach problem approximately from a given
% initial state using a host of techniques
% =============================================================================
%
% SReachPoint computes an approximation to the first/terminal stochastic reach
% problems. 
%
% This function can (approximately) solve two stochastic reachability
% problems:
%
% 1. First hitting-time stochastic reachability problem:
%
%     maximize Prob( \cup_{i=1}^N {\cap_{t=0}^{i-1} 
%                               x_t lies in Safe_t\TargetHyp, x_i\in\TargetHyp})
%     subject to
%           dynamics and bounds on control
%
% 2. Terminal hitting-time stochastic reachability problem (stochastic
% reachability of a target tube):
%
%     maximize Prob( \cap_{i=1}^N x_t lies in Safe_t)
%     subject to
%           dynamics and bounds on control
%
% Using the theory discussed in,
% 
% A. P. Vinod and M. Oishi, HSCC 2019 TODO
%
% We can underapproximate the first hitting-time problem by computing a
% finite maximum of a series of stochastic reachability of a target tube with 
% varying time-horizon. 
%
% For the affine controller synthesis problem, we relax hard bounds on the
% control to a user-specified bound on the probability that the affine
% controller
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
%    SReachTool function: SReachPointCcO
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
%    SReachTool function: SReachPointCcA
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
%    SReachTool function: SReachPointGpO
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
%    SReachTool function: SReachPointScO TODO
%    Paper              : Lesser, Oishi, Erwin TODO.
%
%
% USAGE: TODO
%
% =============================================================================
%
% [approx_reach_prob, opt_controller, varargout] = SReachPoint(prob_str,...
%    method_str, sys, initial_state, safety_tube, options)
% 
% Inputs:
% -------
%   prob_str     - String specifying the problem of interest. For each case, we
%                  compute the optimal value function that maps initial states
%                  to different maximal reach probabilities
%                      1. 'first' : Stay within the safety_tube and reach the
%                                   target set early if possible
%                      2. 'term' : Stay within the safety_tube
%   method_str   - Solution technique to be used.
%                      'chance-open'  -- Convex chance-constrained approach for
%                                        an open-loop controller synthesis
%                      'chance-affine'-- Convex chance-constrained approach for
%                                        an affine controller synthesis
%                      'genzps-open'  -- Genz's algorithm + Patternsearch
%                      'scenario-open'-- Scenario-based 
%   sys          - System description (LtvSystem/LtiSystem object)
%   initial_state- Initial state for which the maximal reach probability must be
%                  evaluated (A numeric vector of dimension sys.state_dim)
%   safety_tube  - Collection of (potentially time-varying) safe sets that
%                  define the safe states (TargetTube object)
%   options      - Collection of user-specified options for each of the solution
%                  (Matlab struct created using SReachPointOptions)
%
% Outputs:
% --------
%   approx_reach_prob 
%               - Approximation (underapproximation, in some cases) to the
%                 first/terminal stochastic reach problem
%   opt_input_vec, 
%     opt_input_gain
%               - Controller U=MW+d for a concatenated input vector 
%                   U = [u_0; u_1; ...; u_{N-1}] and concatenated disturbance
%                   vector W=[w_0; w_1; ...; w_{N-1}]. 
%                   - opt_input_gain: Affine controller gain matrix of dimension
%                       (sys.input_dim*N) x (sys.dist_dim*N)
%                   - opt_input_vec: Open-loop controller: column vector dimension
%                       (sys.input_dim*N) x 1
%                 The feedback gain matrix M is set to [] for methods that look
%                 for open-loop controllers. 
%   risk_alloc_state 
%               - [Available only for 'chance-X'] Risk allocation for the
%                 state constraints
%   risk_alloc_input
%               - [Available only for 'chance-affine'] Risk allocation for the
%                 input constraints
%
% Notes:
% * SReachPoint() will call this function internally using the default
%     values if SReachPointOptions()-based options is not explicitly provided
%     to SReachPoint().
% * See @LtiSystem/getConcatMats for more information about the notation used.
% * If an open_loop policy is desired arranged in increasing time columnwise,
%   use the following command:
%       optimal_open_loop_control_policy = reshape(opt_controller,...
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
    inpar.addRequired('initial_state', @(x) validateattributes(x,...
        {'numeric'}, {'vector'}));
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
        approx_reach_prob = 0;
        time_horizon = length(safety_tube) - 1;
        opt_input_vec = nan(sys.input_dim * time_horizon, 1);
        opt_input_gain = [];        
    elseif strcmpi(prob_str,'term')
        % Ensure that options are provided are appropriate
        options = otherInputHandling(prob_str,method_str, varargin);
        
        % Dependig on method_str call the appropriate solution technique
        switch(lower(method_str))
            case 'genzps-open'
                % Patternsearch and Fourier transform-based open-loop
                % underapproximation of the stochastic reach-avoid problem
                [approx_reach_prob, opt_input_vec] = SReachPointGpO(sys,...
                    initial_state, safety_tube, options);
                 opt_input_gain = [];
            case 'chance-open'
                % Chance-constrained formulation with piecewise-linear 
                % approximations to compute open-loop controller (LP)
                [approx_reach_prob, opt_input_vec, risk_alloc_state] =...
                    SReachPointCcO(sys, initial_state, safety_tube, options);
                 opt_input_gain = [];
                 varargout{1} = risk_alloc_state;
            case 'scenario-open'
                % Chance-constrained formulation with piecewise-linear 
                % approximations to compute open-loop controller (LP)
                [approx_reach_prob, opt_input_vec] = SReachPointCcO(sys,...
                    initial_state, safety_tube, options);
                 opt_input_gain = [];
            case 'chance-affine'
                % Chance-constrained formulation with piecewise-linear 
                % approximations to compute affine-loop controller (SOC program)
                [approx_reach_prob_affine, opt_input_vec, opt_input_gain,...
                    risk_alloc_state, risk_alloc_input] = SReachPointCcA(sys,...
                        initial_state, safety_tube, options);
                approx_reach_prob = approx_reach_prob_affine; %TODO: Do transform
                varargout{1} = risk_alloc_state;
                varargout{2} = risk_alloc_input;
        end
        if approx_reach_prob<0
            approx_reach_prob = 0;
            opt_input_vec = opt_input_vec;
            opt_input_gain = [];
            if strcmpi(method_str,'chance-open') ||...
                    strcmpi(method_str,'chance-affine')    %TODO: Do MILP
                warn_str = ['Note that ''chance-open/chance-affine'' works ',...
                    'only if the maximal reach probability is above 0.5.'];
            else
                warn_str = 'The problem may be infeasible';
            end
            % Need to add sprintf (despite MATLAB's editor warning) for new line
            warning('SReachTools:runtime',sprintf(['%s returned a ',...
                'trivial lower bound.\n%s'], method_str, warn_str));             
        end
    elseif strcmpi(prob_str,'first')
        % Ensure that options, target_hyperplane are provided are appropriate
        [options, target_hyperplane] =otherInputHandling(prob_str,method_str,... 
            varargin, sys.state_dim);
        
        % No need to do any control, the probability is strictly one
        if target_hyperplane.contains(initial_state)
            approx_reach_prob = 1;
            opt_input_vec = nan(sys.input_dim * time_horizon, 1);
            opt_input_gain = [];
        else
            approx_reach_prob_vec = zeros(1,time_horizon);
            opt_input_vec_cells = cell(1,time_horizon);
            opt_input_gain_cells = cell(1,time_horizon);
            
            % Safety tube intersection with the hyperplane:
            % Complement of ax <= b is ax > b is -ax <-b
            target_hyperplane_complement = Polyhedron('H',-target_hyperplane.H);
            safety_tube_minus_target_set = TargetTube('intersect',...
                safety_tube, target_hyperplane_complement);

            % Compute the series of terminal problems (stochastic
            % reachability of a target tube of appropriate lengths)
            prev_warning_state = getSrtWarnings();
            setSrtWarnings('all','off');
            for t=1:time_horizon
                % Collect all the safe set \setminus target_hyperplane into
                % a cell array of polyhedra (TODO: obtain a cell array)
                safe_sets = safety_tube_minus_target_set.extract([1 t]);
                % Construct the reach tube up to t-1 (count in theory starts 
                % from 0) and set the target set at t
                reach_tube_at_t = TargetTube(safe_sets{:}, target_hyperplane);
                options.prob_str = 'term';
                % The terminal subproblem that is to solved for each t
                [approx_reach_prob_vec(t), opt_input_vec_cells{t},...
                    opt_input_gain_cells{t}] = SReachPoint('term',...
                        method_str,sys,initial_state,reach_tube_at_t, options);
            end
            setSrtWarnings('all',prev_warning_state);
            
            % Find the time horizon that does the best
            [approx_reach_prob, opt_time_indx] = max(approx_reach_prob_vec);
            opt_input_vec = opt_input_vec_cells{opt_time_indx};
            opt_input_gain = opt_input_gain_cells{opt_time_indx};            
        end
    else
        % Exhausted all options => prob_str can be first or term only due
        % to input handling
        throw(SrtDevError('Dealing with an unknown problem configuration'));
    end
end

function [options, varargout] = otherInputHandling(prob_str,method_str,...
    additional_args, varargin)
    % input handling for options, [target_hyperplane]
    % for term, all arguments are explicit and output is input santizied options
    % for first, an additional input argument of sys is required, and the
    % output includes an input santizied target_hyperplane

    if strcmpi(prob_str,'term') && length(additional_args) <= 1        
        if isempty(additional_args) || isempty(additional_args{1})
            % No options provided! Create one.
            options = SReachPointOptions(prob_str,method_str);
        else
            % Options provided
            options = additional_args{1};        
        end
    elseif strcmpi(prob_str,'first') && length(additional_args)<= 2
        state_dim = varargin{1};
        if length(additional_args)<= 1
            throwAsCaller(SrtInvalidArgsError(['Expected options (should be',...
                ' left empty for default options) and a target hyperplane']));
        elseif isempty(additional_args{1})
            % No options provided! Create one.
            options = SReachPointOptions(prob_str,method_str);
        else
            options = additional_args{1};
            target_hyperplane = additional_args{2};
            target_hyperplane.minHRep();
        end
        if ~isequal(size(target_hyperplane.H),[1 state_dim+1])
            throwAsCaller(SrtInvalidArgsError(['Inconsistent target ',...
                'hyperplane provided. target_hyperplane.H should be a 1 x ',...
                'sys.state_dim dimensional row vector']));    
        end
        varargout{1} = target_hyperplane;
    else
        % Lots more than needed!
        throwAsCaller(SrtInvalidArgsError('Too many input arguments'));
    end
    %Check if prob_str and method_str are consistent        
    if ~strcmpi(options.prob_str,'term')
        throwAsCaller(...
            SrtInvalidArgsError('Mismatch in prob_str in the options'));
    end
    if ~strcmpi(options.method_str,method_str)
        throwAsCaller(...
            SrtInvalidArgsError('Mismatch in method_str in the options'));
    end        
end
