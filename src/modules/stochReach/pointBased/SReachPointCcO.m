function [lb_stoch_reach, opt_input_vec, risk_alloc_state, varargout] =...
    SReachPointCcO(sys, initial_state, safety_tube, options)
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
%   stoch_reach_prob 
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
    [concat_safety_tube_A, concat_safety_tube_b] = safety_tube.concat([2 time_horizon+1]);

    %% Halfspace-representation of U^N, H, mean_X_sans_input, cov_X_sans_input
    % GUARANTEES: Non-empty input sets (polyhedron)
    [concat_input_space_A, concat_input_space_b] = getConcatInputSpace(sys,...
        time_horizon);
    % GUARANTEES: Compute the input concat and disturb concat transformations
    [~, H, ~] = getConcatMats(sys, time_horizon);
    % GUARANTEES: Gaussian-perturbed LTI system (sys) and well-defined
    % initial_state and time_horizon
    sysnoi = LtvSystem('StateMatrix',sys.state_mat,'DisturbanceMatrix',...
        sys.dist_mat,'Disturbance',sys.dist);
    [mean_X_sans_input, cov_X_sans_input] = SReachFwd('concat-stoch', sysnoi,...
        initial_state, time_horizon);

    
    %% Compute M --- the number of polytopic halfspaces to worry about
    n_lin_state = size(concat_safety_tube_A,1);

    %% Compute \sqrt{h_i^\top * \Sigma_X_no_input * h_i}
    sigma_vector = diag(sqrt(diag(concat_safety_tube_A *...
        cov_X_sans_input * concat_safety_tube_A')));


    %% Obtain the piecewise linear overapproximation of norminvcdf in [0,0.5]
    [invcdf_approx_m, invcdf_approx_c, lb_deltai] =...
        computeNormCdfInvOverApprox(0.5, options.pwa_accuracy, n_lin_state);

    %% Solve the feasibility problem
    cvx_begin quiet
        variable U_vector(sys.input_dim * time_horizon,1);
        variable mean_X(sys.state_dim * time_horizon, 1);
        variable deltai(n_lin_state, 1);
        variable norminvover(n_lin_state, 1);
        minimize sum(deltai)
        subject to
            mean_X == mean_X_sans_input + H * U_vector;
            concat_input_space_A * U_vector <= concat_input_space_b;
            for deltai_indx=1:n_lin_state
                norminvover(deltai_indx) >= invcdf_approx_m.*...
                    deltai(deltai_indx) + invcdf_approx_c; 
            end
            concat_safety_tube_A * mean_X + sigma_vector * norminvover...
                <= concat_safety_tube_b;
            deltai >= lb_deltai;
            deltai <= 0.5;
            sum(deltai) <= 0.5;
    cvx_end

    %% Overwrite the solutions
    if strcmpi(cvx_status, 'Solved')
        lb_stoch_reach = 1-sum(deltai);
        opt_input_vec = U_vector; 
        risk_alloc_state = deltai;
    else
        lb_stoch_reach = -1;
        opt_input_vec = nan(sys.input_dim * time_horizon,1);
        risk_alloc_state = nan(n_lin_state, 1);
    end
    
    %% Create the other info for use if necessary
    extra_info.concat_safety_tube_A = concat_safety_tube_A;
    extra_info.concat_safety_tube_b = concat_safety_tube_b;
    extra_info.concat_input_space_A = concat_input_space_A;
    extra_info.concat_input_space_b = concat_input_space_b;
    extra_info.H = H;
    extra_info.mean_X_sans_input = mean_X_sans_input;
    extra_info.cov_X_sans_input = cov_X_sans_input;
    varargout{1} = extra_info;
end

function otherInputHandling(options)
    if ~(strcmpi(options.prob_str, 'term') &&...
            strcmpi(options.method_str, 'chance-open'))
        throwAsCaller(SrtInvalidArgsError('Invalid options provided'));
    end
end
