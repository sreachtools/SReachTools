function [lb_stoch_reach, opt_input_vec, risk_alloc_state, varargout] =...
    SReachPointCcO(sys, initial_state, safety_tube, desired_accuracy)
% Solve the stochastic reach-avoid problem (lower bound on the probability and
% an open-loop controller synthesis) using chance-constrained convex
% optimization optimization (Internal function --- assumes arguments are all ok)
% =============================================================================
%
% computeCcLowerBoundStochReachAvoidPwlRisk implements the chance-constrained
% convex underapproximation to the terminal hitting-time stochastic reach-avoid
% problem discussed in
%
% K. Lesser, M. Oishi, and R. Erwin, "Stochastic reachability for control of
% spacecraft relative motion," in IEEE Conference on Decision and Control (CDC),
% 2013.
%
% and reformulated in
%
% A. Vinod and M. Oishi, HSCC 2018 TODO
%
% =============================================================================
%   [lb_stoch_reach_avoid, optimal_input_vector] =...
%       computeCcLowerBoundStochReachAvoidPwlRisk( ...
%           sys, ...
%           time_horizon, ...
%           concat_input_space_A, ... 
%           concat_input_space_b, ...
%           concat_safety_tube_A, ... 
%           concat_safety_tube_b, ...
%           H, ...
%           mean_X_sans_input, ...
%           cov_X_sans_input, ...
%           desired_accuracy)
% 
% Inputs:
% -------
%   sys                   - LtiSystem object
%   time_horizon          - Time horizon (N) with the control provided from 0 to N-1
%   concat_input_space_A, 
%    concat_input_space_b - (A,b) Halfspace representation for the
%                            polytope U^{time_horizon} set.        
%   concat_safety_tube_A, 
%    concat_safety_tube_b - (A,b) Halfspace representation for the target tube
%                            from t=1 to time_horizon.  For example, the
%                            terminal reach-avoid problem requires a polytope of
%                            the form safe_set^{time_horizon-1} x safety_set.        
%   H                     - Concatenated input matrix (see
%                            @LtiSystem/getConcatMats for the notation used)
%   mean_X_sans_input     - Mean of X without the influence of the input
%   cov_X_sans_input      - Covariance of X without the influence of the input
%   desired_accuracy      - Desired accuracy for the optimal stochastic
%                           reach-avoid probability
%
% Outputs:
% --------
%   lb_stoch_reach_avoid - Lower bound on the terminal-hitting stochastic
%                          reach avoid problem computed using Fourier
%                          transform and convex optimization
%   optimal_input_vector - Optimal open-loop policy
%                          ((sys.input_dim) *
%                          time_horizon)-dimensional vector 
%                          U = [u_0; u_1; ...; u_N] (column vector)
%
% See also getLowerBoundStochReachAvoid,
% computeCcLowerBoundStochReachAvoidIterRisk, and
% computeFtLowerBoundStochReachAvoid.
%
% Notes:
% ------
% * NOT ACTIVELY TESTED: Builds on other tested functions.
% * MATLAB DEPENDENCY: Uses MATLAB's Statistics and Machine Learning Toolbox.
%                      Needs norminv
% * EXTERNAL DEPENDENCY: Uses MPT3 and CVX (optional)
%                      Needs MPT3 for defining a controlled system and the
%                      definition of the safe and the target (polytopic) sets
%                      Needs CVX to setup a convex optimization problem that
%                      initializes the patternsearch-based optimization. If CVX
%                      is unavailable, the user may provide a guess for the
%                      initialization.
% * See @LtiSystem/getConcatMats for more information about the
%     notation used.
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
    inpar.addRequired('initial_state', @(x) validateattributes(x, {'numeric'},...
        {'vector'}));
    inpar.addRequired('safety_tube',@(x) validateattributes(x,{'TargetTube'},...
        {'nonempty'}));
    inpar.addRequired('desired_accuracy', @(x) validateattributes(x,...
        {'numeric'}, {'scalar','>',0}));

    try
        inpar.parse(sys, initial_state, safety_tube, desired_accuracy);
    catch err
        exc = SrtInvalidArgsError.withFunctionName();
        exc = exc.addCause(err);
        throwAsCaller(exc);
    end

    % Target tubes has polyhedra T_0, T_1, ..., T_{time_horizon}
    time_horizon = length(safety_tube)-1;

    % Get half space representation of the target tube and time horizon
    % skipping the first time step
    [concat_safety_tube_A, concat_safety_tube_b] = safety_tube.concat([2 time_horizon+1]);

    % Construct U^N 
    % GUARANTEES: Non-empty input sets (polyhedron)
    [concat_input_space_A, concat_input_space_b] = getConcatInputSpace(sys, time_horizon);
    % Compute H, mean_X_sans_input, cov_X_sans_input, G for the
    % safety_cost_function definition
    % GUARANTEES: Gaussian-perturbed LTI system (sys) and well-defined
    % init_state and time_horizon
    [H, mean_X_sans_input, cov_X_sans_input] = ...
        getHmatMeanCovForXSansInput(sys, initial_state, time_horizon);

    
    %% Compute M --- the number of polytopic halfspaces to worry about
    n_lin_state = size(concat_safety_tube_A,1);

    %% Compute \sqrt{h_i^\top * \Sigma_X_no_input * h_i}
    sigma_vector = diag(sqrt(diag(concat_safety_tube_A *...
        cov_X_sans_input * concat_safety_tube_A')));


    %% Obtain the piecewise linear overapproximation of norminvcdf in [0,0.5]
    [invcdf_approx_m, invcdf_approx_c, lb_deltai] =...
        computeNormCdfInvOverApprox(0.5, desired_accuracy, n_lin_state);

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