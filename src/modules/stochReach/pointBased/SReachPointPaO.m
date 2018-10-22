function [approx_stoch_reach, opt_input_vec] = SReachPointPaO(sys, ...
    initial_state, safety_tube, options)
% Solve the problem of stochastic reachability of a target tube (a lower bound
% on the maximal reach probability and an open-loop controller synthesis) using
% particle filter control
% =============================================================================
%
% SReachPointPaO implements a mixed-integer linear program-based approximation
% to the stochastic reachability of a target tube problem. This solution is
% based off the particle filter control formulation (for the simpler terminal
% hitting-time stochastic reach-avoid problem) discussed in
%
% K. Lesser, M. Oishi, and R. Erwin, "Stochastic reachability for control of
% spacecraft relative motion," in IEEE Conference on Decision and Control (CDC),
% 2013.
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
%
% =============================================================================
%
% [approx_stoch_reach, opt_input_vec] = SReachPointPaO(sys, initial_state,...
%   safety_tube, options)
%
% Inputs:
% -------
%   sys          - System description (LtvSystem/LtiSystem object)
%   initial_state- Initial state for which the maximal reach probability must be
%                  evaluated (A numeric vector of dimension sys.state_dim)
%   safety_tube  - Collection of (potentially time-varying) safe sets that
%                  define the safe states (Tube object)
%   options      - Collection of user-specified options for 'particle-open'
%                  (Matlab struct created using SReachPointOptions)
%
% Outputs:
% --------
%   approx_stoch_reach 
%               - An approximation of the stochastic reachability of a target
%                   tube problem computed using particle control
%   opt_input_vec
%               - Open-loop controller: column vector of dimension
%                 (sys.input_dim*N) x 1
%
% See also SReachPoint.
%
% Notes:
% * This function requires CVX with Gurobi as the backend solver for optimizing
%   the resulting mixed-integer linear program.
% * See @LtiSystem/getConcatMats for more information about the notation used.
% 
% ============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 
%

    % Input parsing
    inpar = inputParser();
    inpar.addRequired('sys', @(x) validateattributes(x, ...
        {'LtiSystem','LtvSystem'}, {'nonempty'}));
    inpar.addRequired('initial_state', @(x) validateattributes(x, ...
        {'numeric'}, {'vector'}));
    inpar.addRequired('safety_tube',@(x) validateattributes(x, {'Tube'}, ...
        {'nonempty'}));

    try
        inpar.parse(sys, initial_state, safety_tube);
    catch err
        exc = SrtInvalidArgsError.withFunctionName();
        exc = exc.addCause(err);
        throwAsCaller(exc);
    end
    
    % Target tubes has polyhedra T_0, T_1, ..., T_{time_horizon}
    time_horizon = length(safety_tube)-1;

    approx_stoch_reach = -1;
    opt_input_vec = nan(sys.input_dim * time_horizon,1);

    % Requires Gurobi since we are solving a MILP
    [default_solver, solvers_cvx] = cvx_solver;
    if ~(contains(default_solver,'Gurobi') ||...
        any(contains(solvers_cvx,'Gurobi')))
        warning('SReachTools:setup',['SReachPointPaO returns a trivial ', ...
            'result since Gurobi (a MILP solver) was not setup.']);
    else
        % Ensure options is good
        otherInputHandling(options, sys);
        
        % Get half space representation of the target tube and time horizon
        % skipping the first time step
        [concat_safety_tube_A, concat_safety_tube_b] =...
            safety_tube.concat([2 time_horizon+1]);

        % Halfspace-representation of U^N, and matrices Z, H, and G
        % GUARANTEES: Non-empty input sets (polyhedron)
        [concat_input_space_A, concat_input_space_b] = getConcatInputSpace(...
            sys, time_horizon);
        % Compute the input concatenated transformations
        [Z, H, G] = getConcatMats(sys, time_horizon);
        
        % Compute M --- the number of polytopic halfspaces to worry about
        n_lin_state = size(concat_safety_tube_A,1);

        % Compute the stochasticity of the concatenated disturbance random vec
        muW = repmat(sys.dist.parameters.mean,time_horizon,1);
        covW = kron(eye(time_horizon), sys.dist.parameters.covariance);

        % Create realizations of W arranged columnwise
        if options.verbose >= 1
            fprintf('Creating Gaussian random variable realizations....');
        end    
        W_realizations = mvnrnd(muW', covW, options.num_particles)';
        if options.verbose >= 1
            fprintf('Done\n');
        end

        % Implementation of Problem 2 in Lesser CDC 2013
        % Solve the mixed-integer linear program
        if options.verbose >= 1
            fprintf('Setting up CVX problem....');
        end
        cvx_begin
            if options.verbose >= 2
                cvx_quiet false
            else
                cvx_quiet true
            end
            cvx_solver Gurobi
            variable U_vector(sys.input_dim * time_horizon,1);
            variable mean_X(sys.state_dim * time_horizon,options.num_particles);
            variable z(1,options.num_particles) binary;
            maximize sum(z)
            subject to
                mean_X == repmat(Z * initial_state + H * U_vector, ...
                    1, options.num_particles) + G * W_realizations;
                concat_input_space_A * U_vector <= concat_input_space_b;
                concat_safety_tube_A * mean_X <= repmat( ...
                    concat_safety_tube_b, 1, options.num_particles) + ...
                    options.bigM * repmat(1-z,n_lin_state,1);
            if options.verbose >= 1
                fprintf('Done\nParsing and solving the MILP....');
            end
        cvx_end
        if options.verbose >= 1
            fprintf('Done\n');
        end
        %% Overwrite the solutions
        if strcmpi(cvx_status, 'Solved')
            approx_stoch_reach = sum(z)/options.num_particles;
            opt_input_vec = U_vector; 
        end
    end
end

function otherInputHandling(options, sys)
    if ~(strcmpi(options.prob_str, 'term') &&...
            strcmpi(options.method_str, 'particle-open'))
        throwAsCaller(SrtInvalidArgsError('Invalid options provided'));
    end
    if ~strcmpi(sys.dist.type,'Gaussian')
        throwAsCaller(SrtInvalidArgsError(['Expected a Gaussian-perturbed', ...
            'linear system']));
    end
end
