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
%                 tube problem computed using particle control. While it is
%                 expected to lie in [0,1], it is set to -1 in cases where the
%                 CVX optimization fails (cvx_status \not\in {Solved,
%                 Inaccurate/Solved}).
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
%      https://sreachtools.github.io/license/
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
    n_Gurobi_solver = nnz(contains(solvers_cvx,'Gurobi'));
    if n_Gurobi_solver == 0
        warning('SReachTools:runtime',['SReachPointVoO returns a trivial ', ...
            'result since Gurobi (a MILP solver) was not setup.']);
    else
        if ~contains(default_solver, 'Gurobi') && n_Gurobi_solver >= 2
            warning('SReachTools:runtime', sprintf(['SReachPointVoO ', ...
                'requires a MILP solver. Found %d Gurobi solvers.\n', ...
                'Choosing the Gurobi solver bundled with CVX.\nSet the ',...
                'desired Gurobi solver, before calling this function.'], ...
                n_Gurobi_solver));              
        end
        
        otherInputHandling(sys, options);
        
        % Get half space representation of the target tube and time horizon
        % skipping the first time step
        [concat_safety_tube_A, concat_safety_tube_b] =...
            safety_tube.concat([2 time_horizon+1]);
       
        % Normalize the hyperplanes so that Ax <= b + M(1-bin_x) can be 
        % implemented with M = 100 (numerically stable compared to using large
        % M). Here, b<=1e-3 or 1
        [concat_safety_tube_A, concat_safety_tube_b] =...
            normalizeForParticleControl(concat_safety_tube_A,...
                concat_safety_tube_b);

        % Halfspace-representation of U^N, and matrices Z, H, and G
        % GUARANTEES: Non-empty input sets (polyhedron)
        [concat_input_space_A, concat_input_space_b] = getConcatInputSpace(...
            sys, time_horizon);
        % Compute the input concatenated transformations
        [Z, H, G] = getConcatMats(sys, time_horizon);
        
        % Compute M --- the number of polytopic halfspaces to worry about
        n_lin_state = size(concat_safety_tube_A,1);

        if options.verbose >= 1
            fprintf('Required number of particles: %d\n', ... 
                options.n_particles);
            fprintf('Creating random variable realizations....');
        end        
        % Compute the stochasticity of the concatenated disturbance random vec
        W = concat(sys.dist, time_horizon);        
        % Create realizations of W arranged columnwise
        W_realizations = W.getRealizations(options.n_particles);
        if options.verbose >= 1
            fprintf('Done\n');
        end
        
        % Implementation of Problem 2 in Lesser CDC 2013
        % Solve the mixed-integer linear program
        if options.verbose >= 1
            if options.verbose >= 2
                fprintf('Objective value needs to be scaled by %1.3f\n',...
                    1/options.n_particles);
            end
            fprintf('Setting up CVX problem....');
        end
        cvx_begin
            if options.verbose >= 2
                cvx_quiet false
            else
                cvx_quiet true
            end
            if ~contains(default_solver,'Gurobi')
                cvx_solver Gurobi;
            end
            variable U_vector(sys.input_dim * time_horizon,1);
            variable mean_X(sys.state_dim * time_horizon,options.n_particles);
            variable bin_x(1,options.n_particles) binary;            
            maximize sum(bin_x)/options.n_particles
            subject to
                mean_X == repmat(Z * initial_state + H * U_vector, ...
                    1, options.n_particles) + G * W_realizations;
                concat_input_space_A * U_vector <= concat_input_space_b;
                concat_safety_tube_A * mean_X <= repmat( ...
                    concat_safety_tube_b, 1, options.n_particles) + ...
                    options.bigM * repmat(1-bin_x,n_lin_state,1);
            if options.verbose >= 1
                fprintf('Done\nParsing and solving the MILP....');
            end
        cvx_end
        if options.verbose >= 1
            fprintf('Done\n');
        end
        %% Overwrite the solutions
        switch cvx_status
            case {'Solved', 'Inaccurate/Solved'}
                approx_stoch_reach = sum(bin_x)/options.n_particles;
                opt_input_vec = U_vector; 
            otherwise
                
        end
    end
end

function otherInputHandling(sys, options)
    % 1. Get the correct options
    % 2. Check if the system is stochastic
    if ~(strcmpi(options.prob_str, 'term') &&...
            strcmpi(options.method_str, 'particle-open'))
        throwAsCaller(SrtInvalidArgsError('Invalid options provided'));
    end
    if ~isa(sys.dist,'RandomVector')
        throwAsCaller(SrtInvalidArgsError('Expected a stochastic system'));
    end
end
