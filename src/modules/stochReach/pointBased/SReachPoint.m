function [approx_reach_prob, opt_input_vec, opt_input_gain, varargout] =...
    SReachPoint(prob_str, method_str, sys, initial_state, safety_tube, varargin)
% Solve the problem of stochastic reachability of a target tube from a given
% initial state using a host of techniques
% =============================================================================
%
% SReachPoint can (approximately) solve the problem of stochastic reachability 
% of a target tube from a given initial state,
%
%     maximize Prob( \cap_{i=1}^N x_t lies in Safe_t)
%     subject to
%           dynamics and bounds on control
%
% This function can compute an underapproximation to the optimal value of the
% above problem as well as synthesize controllers:
% 1. open-loop controller that satisfies the prespecified hard control bounds,
% 2. affine disturbance feedback controller that satisfies the hard control
%       bounds with a likelihood above the user-specified threshold.
% This function is a compilation of various techniques proposed in the
% literature.
%
% The problem of stochastic reachability of a target tube is discussed in 
%
% A. Vinod and M. Oishi, "Stochastic reachability of a target tube: Theory and
% computation," IEEE Transactions in Automatic Control, 2018 (submitted) URL:
% https://arxiv.org/pdf/1810.05217.pdf.
%
% It subsumes the terminal hitting-time stochastic reach-avoid problem (Summers
% and Lygeros, Automatica, 2010) as well as the stochastic viability problem
% (Abate et. al, Automatica, 2008).
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
%    Paper              : a. K. Lesser, M. Oishi, and R. Erwin, "Stochastic
%                            reachability for control of spacecraft relative
%                            motion," In Proc. IEEE Conf. Dec. & Ctrl., 2013.
%                         b. A. Vinod and M. Oishi. Affine controller synthesis
%                            for stochastic reachability via difference of
%                            convex programming. In Proc. Conf. of Dec. & Ctrl.,
%                            2019. (submitted).
%                            https://hscl.unm.edu/affinecontrollersynthesis/
%
% 2. Difference-of-convex chance-constrained-based approach (chance-affine):
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
%                         underapproximation problem due to non-convexity
%                         established by the difference of convex formulation
%    Approximation      : Guaranteed underapproximation
%    SReachTool function: SReachPointCcA
%    Dependency (EXT)   : CVX
%    Paper              : A. Vinod and M. Oishi. Affine controller synthesis for
%                         stochastic reachability via difference of convex
%                         programming. In Proc. Conf. of Dec. & Ctrl., 2019.
%                         (submitted).
%                         https://hscl.unm.edu/affinecontrollersynthesis/
%
% 3. Uniform risk allocation-based approach (chance-affine-uni):
%
%    High-level desc.   : Use Boole's inequality, Gaussian random vector,
%                         hyperbolic constraints-to-second order cone constraint
%                         reformulation. By decoupling the risk allocation
%                         from controller synthesis, this approach solves a
%                         series of SOCP problems guided by a bisection on
%                         the uniform risk allocation (See Notes)
%    Controller type    : A history-dependent affine controller that satisfies
%                         softened input constraints (controller satisfies the
%                         hard input bounds upto a user-specified probabilistic
%                         threshold)
%    Optimality         : Suboptimal affine controller for the
%                         underapproximation problem
%    Approximation      : Guaranteed underapproximation
%    SReachTool function: SReachPointCcAu
%    Dependency (EXT)   : CVX
%    Paper              : 1) M. Vitus and C. Tomlin, "On feedback design and 
%                            risk allocation in chance constrained control", In 
%                            Proc. Conf. Dec. & Ctrl., 2011.
%                         2) A. Vinod and M. Oishi. Affine controller synthesis
%                            for stochastic reachability via difference of 
%                            convex programming. In Proc. Conf. of Dec. & Ctrl., 
%                            2019 (submitted).
%                            https://hscl.unm.edu/affinecontrollersynthesis/
%
% 4. Fourier transform + Patternsearch (genzps-open):
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
% 5. Particle control approach (particle-open):
%
%    High-level desc.   : Sample particles based on the additive noise and solve
%                         a mixed-integer linear program to make the maximum
%                         number of particles satisfy the reachability objective
%    Approximation      : No direct approximation guarantees. Accuracy improves
%                         as the number of particles considered increases.
%    Controller type    : Open-loop controller that satisfies the hard input
%                         bounds
%    Optimality         : Optimal (w.r.t particles drawn) open-loop controller
%                         for the underapproximation problem 
%    Dependency (EXT)   : CVX, Gurobi
%    SReachTool function: SReachPointPaO
%    Paper              : K. Lesser, M. Oishi, and R. Erwin, "Stochastic
%                         reachability for control of spacecraft relative
%                         motion," In Proc. IEEE Conf. Dec. & Ctrl., 2013.
%
% 6. Particle control-based approach with undersampling via Voronoi partitions 
%    (voronoi-open):
%
%    High-level desc.   : Sample particles based on the additive noise and solve
%                         a mixed-integer linear program to make the maximum
%                         number of particles satisfy the reachability objective
%                         In addition, we use Voronoi partition to
%                         drastically improve the tractability while
%                         preserving the underapproximation quality
%    Approximation      : Overapproximation bounded above (in probability) by a
%                         user-specified tolerance
%    Controller type    : Open-loop controller that satisfies the hard input
%                         bounds
%    Optimality         : Optimal (w.r.t particles drawn) open-loop controller
%                         for the underapproximation problem 
%    Dependency (EXT)   : CVX, Gurobi
%    SReachTool function: SReachPointVoO
%    Paper              : H. Sartipizadeh, A. Vinod, B. Acikmese, and M. Oishi, 
%                         "Voronoi Partition-based Scenario Reduction for Fast 
%                         Sampling-based Stochastic Reachability Computation of
%                         LTI Systems", In Proc. Amer. Ctrl. Conf., 2019
%
% See also examples/cwhSReachPointDemo.m and examples/dubinsSReachPointDemo.m.
%
% =============================================================================
%
% [approx_reach_prob, opt_input_vec, opt_input_gain] = SReachPoint(prob_str,...
%   method_str, sys, initial_state, safety_tube, [options])
% 
% Inputs:
% -------
%   prob_str     - String specifying the problem of interest. For each case, we
%                  compute the optimal value function that maps initial states
%                  to different maximal reach probabilities
%                      1. 'term' : Stay within the safety_tube
%   method_str   - Solution technique to be used.
%                      'chance-open'  -- Convex chance-constrained approach for
%                                        an open-loop controller synthesis;
%                                        performs risk allocation and
%                                        controller synthesis
%                                        simultaneously
%                      'chance-affine'-- Difference-of-convex chance-constrained 
%                                        approach for an affine controller 
%                                        synthesis; Allocates risk and 
%                                        synthesizes controller simultaneously
%                      'chance-affine-uni'
%                                     -- Uniform risk allocation approach for an 
%                                        affine controller synthesis; Decouples
%                                        the risk allocation problem from
%                                        controller synthesis
%                      'genzps-open'  -- Genz's algorithm + Patternsearch for an
%                                        open-loop controller synthesis
%                      'particle-open'-- Particle control-based approach for an
%                                        open-loop controller synthesis
%                      'voronoi-open' -- Voronoi undersampling of Particle
%                                        control-based approach for an open-loop
%                                        controller synthesis
%   sys          - System description (LtvSystem/LtiSystem object)
%   initial_state- Initial state for which the maximal reach probability must be
%                  evaluated (A numeric vector of dimension sys.state_dim)
%   safety_tube  - Collection of (potentially time-varying) safe sets that
%                  define the safe states (Tube object)
%   options      - [Optional] Collection of user-specified options for each of
%                  the solution (Matlab struct created using SReachPointOptions)
%
% Outputs:
% --------
%   approx_reach_prob 
%               - Approximation (underapproximation, in some cases) to the
%                 stochastic reach problem
%   opt_input_vec, 
%     opt_input_gain
%               - Controller U=MW+d for a concatenated input vector 
%                   U = [u_0; u_1; ...; u_{N-1}] and concatenated disturbance
%                   vector W=[w_0; w_1; ...; w_{N-1}]. 
%                   - opt_input_gain: Affine controller gain matrix of dimension
%                       (sys.input_dim*N) x (sys.dist_dim*N)
%                   - opt_input_vec: Open-loop controller: column vector 
%                     dimension
%                       (sys.input_dim*N) x 1
%                 The feedback gain matrix M is set to [] for methods that look
%                 for open-loop controllers. 
%   risk_alloc_state 
%               - [Available only for 'chance-X'] Risk allocation for the
%                 state constraints
%   risk_alloc_input
%               - [Available only for 'chance-affine(-X)'] Risk allocation for 
%                 the input constraints
%   extra_info  - [Available only for 'voronoi-X'] MATLAB struct
%                 containing additional info about the Voronoi partition-based 
%                 undersampling used for tractable particle control approach
%
% Notes:
% * SReachPoint() will call SReachPointOptions() internally if
%       SReachPointOptions()-based options is not explicitly provided to
%       SReachPoint(). This will set the algorithm to default options.
% * 'chance-affine' and 'chance-affine-uni' requires an explicit declaration of 
%   the options from SReachPointOptions() to specify the threshold on the 
%   chance-constraint relaxation of the input bounds.
% * 'chance-affine' provides the best lower bound at the cost of additional
%   compute time, when compared to 'chance-affine-uni'
% * See @LtiSystem/getConcatMats for more information about the notation used.
% * If an open_loop policy is desired arranged in increasing time columnwise,
%   use the following command:
%       optimal_open_loop_control_policy = reshape(opt_controller, ...
%           sys.input_dim, time_horizon);
% 
% =============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://sreachtools.github.io/license/
%
%

    % Input parsing
    valid_prob = {'term'}; % TODO-first: 'first',
    valid_method= {'chance-open', 'chance-affine', 'genzps-open',...
        'particle-open', 'voronoi-open', 'chance-affine-uni'};

    inpar = inputParser();
    inpar.addRequired('prob_str', @(x) any(validatestring(x,valid_prob)));
    inpar.addRequired('method_str', @(x) any(validatestring(x,valid_method)));
    inpar.addRequired('sys', @(x) validateattributes(x, ...
        {'LtiSystem','LtvSystem'}, {'nonempty'}));
    inpar.addRequired('initial_state', @(x) validateattributes(x, ...
        {'numeric'}, {'vector'}));
    inpar.addRequired('safety_tube',@(x) validateattributes(x,{'Tube'}, ...
        {'nonempty'}));

    try
        inpar.parse(prob_str, method_str, sys, initial_state, safety_tube);
    catch err
        exc = SrtInvalidArgsError.withFunctionName();
        exc = exc.addCause(err);
        throwAsCaller(exc);
    end

    if ~isa(sys.dist,'RandomVector')
        throwAsCaller(SrtInvalidArgsError('Expected a stochastic system'));
    end
        
    time_horizon = length(safety_tube) - 1;
        
    % Check if safe set contains the initial state
    if ~safety_tube(1).contains(initial_state)
        % Stochastic reach-avoid probability is zero and no admissible open-loop
        % policy exists, if given an unsafe initial state
        approx_reach_prob = 0;
        opt_input_vec = nan(sys.input_dim * time_horizon, 1);
        opt_input_gain = [];        
    elseif strcmpi(prob_str,'term')
        % 1. Ensure that options are provided are appropriate
        % 2. Ensure that prob_str, method_str are consistent with options
        %    provided, if any
        options = otherInputHandling(prob_str,method_str, varargin);
        
        % Dependig on method_str call the appropriate solution technique
        switch(lower(method_str))
            case 'genzps-open'
                % Patternsearch and Fourier transform-based open-loop
                % underapproximation of the stochastic reach-avoid problem
                [approx_reach_prob, opt_input_vec] = SReachPointGpO(sys, ...
                    initial_state, safety_tube, options);
                 opt_input_gain = [];
            case 'chance-open'
                % Chance-constrained formulation with piecewise-linear 
                % approximations to compute open-loop controller (LP)
                if nargout <= 4 % Four because opt_input_gain = []
                    [approx_reach_prob, opt_input_vec, risk_alloc_state] = ...
                        SReachPointCcO(sys, initial_state, safety_tube,options);
                    opt_input_gain = [];
                    varargout{1} = risk_alloc_state;
                elseif nargout == 5
                    [approx_reach_prob, opt_input_vec, risk_alloc_state, ...
                        extra_info] = SReachPointCcO(sys, initial_state, ...
                            safety_tube, options);
                    opt_input_gain = [];
                    varargout{1} = risk_alloc_state;
                    varargout{2} = extra_info;
                else
                    throw(SrtRuntimeError('Too many output arguments'));
                end
            case 'particle-open'
                % Particle control-based approach to compute open-loop
                % controller (MILP)
                [approx_reach_prob, opt_input_vec] = SReachPointPaO(sys, ...
                    initial_state, safety_tube, options);
                 opt_input_gain = [];
            case 'voronoi-open'
                % Undersampled particle control approach to compute
                % open-loop controller (MILP)
                if nargout <= 3 % Three because opt_input_gain = []
                    [approx_reach_prob, opt_input_vec] = SReachPointVoO(sys, ...
                        initial_state, safety_tube, options);
                    opt_input_gain = [];
                elseif nargout == 4
                    [approx_reach_prob, opt_input_vec, extra_info] =...
                        SReachPointVoO(sys, initial_state, safety_tube, ...
                        options);
                     opt_input_gain = [];
                     varargout{1} = extra_info;
                else
                    throw(SrtRuntimeError('Too many output arguments'));
                end                
            case 'chance-affine'
                % Chance-constrained formulation with piecewise-linear 
                % approximations to compute affine-loop controller (SOC program)
                if nargout <= 5
                    [approx_reach_prob, opt_input_vec, opt_input_gain, ...
                        risk_alloc_state, risk_alloc_input] = SReachPointCcA(...
                            sys, initial_state, safety_tube, options);
                    varargout{1} = risk_alloc_state;
                    varargout{2} = risk_alloc_input;
                else
                    throw(SrtRuntimeError('Too many output arguments'));
                end                
            case 'chance-affine-uni'
                % Chance-constrained formulation with piecewise-linear 
                % approximations to compute affine-loop controller (SOC program)
                if nargout <= 5
                    [approx_reach_prob, opt_input_vec, opt_input_gain, ...
                        risk_alloc_state, risk_alloc_input] = ...
                            SReachPointCcAu(sys, initial_state, safety_tube, ...
                                options);
                    varargout{1} = risk_alloc_state;
                    varargout{2} = risk_alloc_input;
                else
                    throw(SrtRuntimeError('Too many output arguments'));
                end                
            otherwise
                throw(SrtInternalError('Internal function not setup!'));
        end
        if approx_reach_prob<0
            approx_reach_prob = 0;
            opt_input_vec = nan(sys.input_dim * time_horizon, 1);
            opt_input_gain = [];
            if strcmpi(method_str,'chance-open') ||...
                    strcmpi(method_str,'chance-affine')  %TODO: Do MILP
                warn_str = sprintf(['Note that ''%s''', ...
                    'works only if all the individual hyperplanes is ', ...
                    'satisfied with a probability above 0.5.'], method_str);
            elseif strcmpi(method_str,'chance-affine-uni')
                warn_str = ['Uniform risk allocation is not possible. Try ', ...
                    'using ''chance-affine''.'];
            else
                warn_str = 'The problem may be infeasible';
            end
            % Need to add sprintf (despite MATLAB's editor warning) for new line
            warning('SReachTools:runtime',sprintf(['%s returned a ', ...
                'trivial lower bound.\n%s'], method_str, warn_str));             
        end
    elseif strcmpi(prob_str,'first')
        % Ensure that options, target_hyperplane are provided are appropriate
        [options, target_hyperplane] = otherInputHandling(prob_str, ...
            method_str, varargin, sys.state_dim);
        
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
            safety_tube_minus_target_set = Tube('intersect', ...
                safety_tube, target_hyperplane_complement);

            % Compute the series of terminal problems (stochastic
            % reachability of a target tube of appropriate lengths)
            prev_warning_state = getSrtWarnings();
            setSrtWarnings('all','off');
            for t = 1:time_horizon
                % Collect all the safe set \setminus target_hyperplane into
                % a cell array of polyhedra (TODO: obtain a cell array)
                safe_sets = safety_tube_minus_target_set.extract([1 t]);
                % Construct the reach tube up to t-1 (count in theory starts 
                % from 0) and set the target set at t
                reach_tube_at_t = Tube(safe_sets{:}, target_hyperplane);
                options.prob_str = 'term';
                % The terminal subproblem that is to solved for each t
                [approx_reach_prob_vec(t), opt_input_vec_cells{t}, ...
                    opt_input_gain_cells{t}] = SReachPoint('term', ...
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

function [options, varargout] = otherInputHandling(prob_str,method_str, ...
    additional_args, varargin)
    % 1. Ensure that options are provided are appropriate
    % 2. Ensure that prob_str, method_str are consistent with options
    %    provided, if any
    % input handling for options, [target_hyperplane]
    % for term, all arguments are explicit and output is input santizied options
    % for first, an additional input argument of sys is required, and the
    % output includes an input santizied target_hyperplane

    if strcmpi(prob_str,'term') && length(additional_args) <= 1        
        %% terminal time problem handling
        if isempty(additional_args) || isempty(additional_args{1})
            if strcmpi(method_str, 'chance-affine') || ...
                    strcmpi(method_str, 'chance-affine-uni') 
                throwAsCaller(SrtInvalidArgsError(sprintf(['Default options',...
                    ' for %s requires a user specified parameter.\n', ...
                    'Expected max_input_viol_prob, the maximum allowed ', ...
                    'likelihood of violating the input constraints.'], ...
                    method_str)));
            else
                % No options provided! Create one.
                options = SReachPointOptions(prob_str,method_str);
            end
        else
            % Options provided
            options = additional_args{1};        
        end
    elseif strcmpi(prob_str,'first') && length(additional_args)<= 2
        %% First time problem handling
        state_dim = varargin{1};
        if length(additional_args)<= 1
            throwAsCaller(SrtInvalidArgsError(['Expected options (should ', ...
                'be left empty for default options) and a target hyperplane']));
        elseif isempty(additional_args{1})
            % No options provided! Create one.
            options = SReachPointOptions(prob_str,method_str);
        else
            options = additional_args{1};
            target_hyperplane = additional_args{2};
            target_hyperplane.minHRep();
        end
        if ~isequal(size(target_hyperplane.H),[1 state_dim+1])
            throwAsCaller(SrtInvalidArgsError(['Inconsistent target ', ...
                'hyperplane provided. target_hyperplane.H should be a 1 x ', ...
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
