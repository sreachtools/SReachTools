function [lb_stoch_reach, opt_input_vec] = SReachPointPaO(sys, initial_state,...
    safety_tube, options)
% Solve the stochastic reach-avoid problem (lower bound on the probability and
% open-loop controller synthesis) using chance-constrained convex optimization
% =============================================================================
%
% SReachPointCcO implements a chance-constrained convex underapproximation to
% the stochastic reachability of a target tube problem. This solution is based
% off the formulation (A) discussed in
%
% K. Lesser, M. Oishi, and R. Erwin, "Stochastic reachability for control of
% spacecraft relative motion," in IEEE Conference on Decision and Control (CDC),
% 2013.
%
% This function implements a convex solver-friendly piecewise-affine restriction
% of the formulation (A), as discussed in
%
% A. Vinod and M. Oishi, HSCC 2018 TODO
%
% =============================================================================
%
% [lb_stoch_reach, opt_input_vec, risk_alloc_state, varargout] =...
%    SReachPointCcO(sys, initial_state, safety_tube, options)
%
% Inputs:
% -------
%   sys          - System description (LtvSystem/LtiSystem object)
%   initial_state- Initial state for which the maximal reach probability must be
%                  evaluated (A numeric vector of dimension sys.state_dim)
%   safety_tube  - Collection of (potentially time-varying) safe sets that
%                  define the safe states (TargetTube object)
%   options      - Collection of user-specified options for 'chance-open'
%                  (Matlab struct created using SReachPointOptions)
%
% Outputs:
% --------
%   lb_stoch_reach 
%               - Lower bound on the stochastic reachability of a target tube
%                   problem computed using convex chance constraints and
%                   piecewise affine approximation
%   opt_input_vec
%               - Open-loop controller: column vector of dimension
%                 (sys.input_dim*N) x 1
%   risk_alloc_state 
%               - Risk allocation for the state constraints
%   extra_info  - [Optional] Useful information to construct the
%                   reachability problem | Used by 'genzps-open' to avoid 
%                   unnecessary recomputation
%                 Matlab struct with members --- concat_safety_tube_A,
%                   concat_safety_tube_b, concat_input_space_A,
%                   concat_input_space_b, H, mean_X_sans_input,
%                   cov_X_sans_input;
%
% See also SReachPoint.
%
% Notes:
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
    inpar.addRequired('sys', @(x) validateattributes(x,...
        {'LtiSystem','LtvSystem'}, {'nonempty'}));
    inpar.addRequired('initial_state', @(x) validateattributes(x, ...
        {'numeric'}, {'vector'}));
    inpar.addRequired('safety_tube',@(x) validateattributes(x,{'TargetTube'},...
        {'nonempty'}));

    try
        inpar.parse(sys, initial_state, safety_tube);
    catch err
        exc = SrtInvalidArgsError.withFunctionName();
        exc = exc.addCause(err);
        throwAsCaller(exc);
    end

    % Ensure options is good
    otherInputHandling(options);
    
    % Target tubes has polyhedra T_0, T_1, ..., T_{time_horizon}
    time_horizon = length(safety_tube)-1;

    % Get half space representation of the target tube and time horizon
    % skipping the first time step
    [concat_safety_tube_A, concat_safety_tube_b] =...
        safety_tube.concat([2 time_horizon+1]);

    %% Halfspace-representation of U^N, H, mean_X_sans_input, cov_X_sans_input
    % GUARANTEES: Non-empty input sets (polyhedron)
    [concat_input_space_A, concat_input_space_b] = getConcatInputSpace(sys,...
        time_horizon);
    % Compute the input concatenated transformations
    [Z, H, G] = getConcatMats(sys, time_horizon);
    
    %% Compute M --- the number of polytopic halfspaces to worry about
    n_lin_state = size(concat_safety_tube_A,1);

    if strcmpi(sys.dist.type,'Gaussian')
        muW = repmat(sys.dist.parameters.mean,time_horizon,1);
        covW = kron(eye(time_horizon), sys.dist.parameters.covariance);
        if options.verbose >= 1
            fprintf('Creating Gaussian random variable realizations....');
        end    
        % Create realizations of W arranged columnwise
        W_realizations = mvnrnd(muW', covW, options.num_particles)';
        if options.verbose >= 1
            fprintf('Done\n');
        end
    else
        throw(SrtInvalidArgsError('Expected a Gaussian-perturbed linear system'));
    end
    %% Solve the feasibility problem
    if options.verbose >= 1
        fprintf('Setting up CVX problem....');
    end
    cvx_begin
        if options.verbose >= 2
            cvx_quiet false
        else
            cvx_quiet true
        end
        variable U_vector(sys.input_dim * time_horizon,1);
        variable mean_X(sys.state_dim * time_horizon,options.num_particles);
        variable z(1,options.num_particles) binary;
        minimize sum(z)
        subject to
            mean_X == repmat(Z * initial_state + H * U_vector,...
                1, options.num_particles) + G * W_realizations;
            concat_input_space_A * U_vector <= concat_input_space_b;
            concat_safety_tube_A * mean_X <= repmat(concat_safety_tube_b, 1,...
                options.num_particles) + options.bigM * repmat(z,n_lin_state,1);
        if options.verbose >= 1
            fprintf('Done\nParsing and solving the MILP....');
        end
    cvx_end
    if options.verbose >= 1
        fprintf('Done\n');
    end
    %% Overwrite the solutions
    if strcmpi(cvx_status, 'Solved')
        lb_stoch_reach = sum(z)/options.num_particles;
        opt_input_vec = U_vector; 
    else
        lb_stoch_reach = -1;
        opt_input_vec = nan(sys.input_dim * time_horizon,1);
    end
end

function otherInputHandling(options)
    if ~(strcmpi(options.prob_str, 'term') &&...
            strcmpi(options.method_str, 'particle-open'))
        throwAsCaller(SrtInvalidArgsError('Invalid options provided'));
    end
end
