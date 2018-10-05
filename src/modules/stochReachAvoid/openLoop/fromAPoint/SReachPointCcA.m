function [stoch_reach_prob, opt_input_vec, opt_input_gain,...
    risk_alloc_state, risk_alloc_input] = SReachPointCcA(sys, initial_state,...
        safety_tube, options)
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
% A. Vinod, V. Sivaramakrishnan, and M. Oishi, CSS-L, 2018 (submitted) TODO
%
% USAGE: This function is intended for internal use as it does not sanitize the
% inputs. Please use getLowerBoundStochReachAvoid instead.
%
% =============================================================================
%   [stoch_reach_prob, opt_input_vec] =...
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
%   stoch_reach_prob - Lower bound on the terminal-hitting stochastic
%                          reach avoid problem computed using Fourier
%                          transform and convex optimization
%   opt_input_vec - Optimal open-loop policy
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
    
    try
        inpar.parse(sys, initial_state, safety_tube);
    catch err
        exc = SrtInvalidArgsError.withFunctionName();
        exc = exc.addCause(err);
        throwAsCaller(exc);
    end

    % Target tubes has polyhedra T_0, T_1, ..., T_{time_horizon}
    time_horizon = length(safety_tube) - 1;

    % Get half space representation of the target tube and time horizon
    % skipping the first time step
    [concat_safety_tube_A, concat_safety_tube_b] = safety_tube.concat(...
        [2 time_horizon+1]);

    % Construct U^N 
    % GUARANTEES: Non-empty input sets (polyhedron)
    [concat_input_space_A, concat_input_space_b] = getConcatInputSpace(sys, time_horizon);
    % GUARANTEES: Compute the input concat and disturb concat transformations
    [~, H, G] = getConcatMats(sys, time_horizon);
    % GUARANTEES: Gaussian-perturbed LTI system (sys) and well-defined
    % init_state and time_horizon
    [H, mean_X_sans_input, ~]=getHmatMeanCovForXSansInput(sys, initial_state,...
        time_horizon);
%     sysnoi = LtvSystem('StateMatrix',sys.state_mat,'DisturbanceMatrix',sys.dist_mat,'Disturbance',sys.dist);
%     [mean_X_sans_input,~] = SReachFwd('concat-stoch', sysnoi, initial_state, time_horizon);
    
    %% Compute M --- the number of polytopic halfspaces to worry about
    n_lin_state = size(concat_safety_tube_A,1);
    n_lin_input = size(concat_input_space_A,1);
    
    %% Covariance of W vector
    cov_concat_disturb = kron(eye(time_horizon),sys.dist.parameters.covariance);
    % Compute a sparse square root of a matrix
    sqrt_cov_concat_disturb = chol(cov_concat_disturb);    
    
    %% Difference of convex penalty approach
    lb_max_state_viol_prob = options.bisect_lb;
    ub_max_state_viol_prob = options.bisect_ub;    
    
    while ub_max_state_viol_prob - lb_max_state_viol_prob > options.desired_accuracy
        % Current iteration for viol probability
        max_state_viol_prob=(ub_max_state_viol_prob + lb_max_state_viol_prob)/2;
        if options.verbose >= 1
            fprintf(['Delta_x = %1.3f | Gap for bisection = %1.3f |'...
                     ' Max gap = %1.3f\n\n'],...
                max_state_viol_prob,...
                ub_max_state_viol_prob - lb_max_state_viol_prob,...
                options.desired_accuracy);
        end
        % Obtain the piecewise linear overapproximation of norminvcdf
        % TODO: Translate desired_accuracy in delta to accuracy in terms of
        % phiinv
        [invcdf_approx_m, invcdf_approx_c, lb_deltai] =...
            computeNormCdfInvOverApprox(max_state_viol_prob,...
                options.desired_accuracy,...
                n_lin_state);
        
        [sol_struct_temp, dc_feas] = differenceOfConvexApproach(...
            max_state_viol_prob,...
            options.max_input_viol_prob,... 
            sys,... 
            time_horizon,...
            n_lin_state,...
            n_lin_input,... 
            mean_X_sans_input,...
            H,...
            G,...
            concat_safety_tube_A,...
            concat_safety_tube_b,...
            concat_input_space_A,...
            concat_input_space_b,...
            sqrt_cov_concat_disturb,...
            invcdf_approx_m,...
            invcdf_approx_c,...
            lb_deltai,...
            options);
        if dc_feas
            % Decrease the upper bound
            ub_max_state_viol_prob = max_state_viol_prob;
            % Update the sol_struct
            sol_struct = sol_struct_temp;
            if options.verbose >= 1
                disp('Feasible solution found! Decreasing (tightening) Delta_x');
            end
        elseif abs(ub_max_state_viol_prob - 1) < eps
            % Infeasible for 0.5 itself
            % Exit from the while loop by setting the gap to be zero
            ub_max_state_viol_prob = lb_max_state_viol_prob;
            % Tell SReachPoint that no solution was found
            sol_struct.stoch_reach_prob = -1;
            sol_struct.opt_input_vec = nan(sys.input_dim * time_horizon,1);
            sol_struct.opt_input_gain = [];
            sol_struct.risk_alloc_state = nan(n_lin_state,1);
            sol_struct.risk_alloc_input = nan(n_lin_input,1);
            if options.verbose >= 1
                disp('Infeasible for Delta_x = 0.5! No solution expected.');
            end            
        else
            % Increase the lower bound
            lb_max_state_viol_prob = max_state_viol_prob;
            if options.verbose >= 1
                disp('Infeasible! Increasing (relaxing) Delta_x');
            end
        end
    end
    
    % Unpackage sol_struct
    stoch_reach_prob = sol_struct.stoch_reach_prob;
    opt_input_vec = sol_struct.opt_input_vec;
    opt_input_gain = sol_struct.opt_input_gain;
    risk_alloc_state = sol_struct.risk_alloc_state;
    risk_alloc_input = sol_struct.risk_alloc_input;
end

%% Difference of convex subproblem implemented in CVX
function [sol_struct, dc_feas] = differenceOfConvexApproach(...
            max_state_viol_prob,...
            max_input_viol_prob,... 
            sys,... 
            time_horizon,...
            n_lin_state,...
            n_lin_input,... 
            mean_X_sans_input,...
            H,...
            G,...
            concat_safety_tube_A,...
            concat_safety_tube_b,...
            concat_input_space_A,...
            concat_input_space_b,...
            sqrt_cov_concat_disturb,...
            invcdf_approx_m,...
            invcdf_approx_c,...
            lb_deltai,...
            options)
    
    % Counter for the iterations        
    iter_count = 0;      
    
    % Initializations for DC iterative algorithm
    obj_curr = Inf;      
    dc_slack_with_tau_curr = Inf;
    
    % First solve 
    n_pwa = length(invcdf_approx_m);
    slack_cc_sqrt_state_iter = zeros(n_lin_state, n_pwa);
    slack_cc_sqrt_input_iter = zeros(n_lin_input, n_pwa);            
    tau_iter = options.tau_initial;

    continue_condition = 1;    
    while continue_condition == 1
        % Store previous iterations
        obj_prev = obj_curr;
        dc_slack_with_tau_prev = dc_slack_with_tau_curr;
            
        % The iteration values are updated at the end of the problem
        cvx_begin quiet
            variable M_matrix(sys.input_dim*time_horizon,sys.dist_dim*time_horizon);
            variable d_vector(sys.input_dim * time_horizon, 1);
            variable mean_X(sys.state_dim * time_horizon, 1);
            % State chance constraint
            variable deltai(n_lin_state, 1) nonnegative;
            variable norm_state_replace_slack(n_lin_state, 1) nonnegative;
            variable slack_cc_sqrt_state(n_lin_state, n_pwa) nonnegative;
            variable slack_reverse_state(n_lin_state, n_pwa) nonnegative;
            % Input chance constraint
            variable gammai(n_lin_input, 1) nonnegative;
            variable norm_input_replace_slack(n_lin_input, 1) nonnegative;
            variable slack_cc_sqrt_input(n_lin_input, n_pwa) nonnegative;
            variable slack_reverse_input(n_lin_input, n_pwa) nonnegative;
            % Minimize slack variable for the norm replacements (epigraph
            % construction) and also the DC programming-based slack constraints
            minimize (sum(norm_input_replace_slack) +...
                        sum(norm_state_replace_slack)...
                            + tau_iter * (sum(sum(slack_reverse_state)) + ...
                                            sum(sum(slack_reverse_input))));
            subject to
                % Causality constraints on M_matrix
                for time_indx = 1:time_horizon - 1
                    M_matrix((time_indx-1)*sys.input_dim + 1:...
                        time_indx*sys.input_dim,...
                        (time_indx-1)*sys.dist_dim+1:end) == 0; 
                end
                % Mean trajectory constraint
                mean_X == mean_X_sans_input + H * d_vector;
                % Risk allocation bounds
                lb_deltai <= deltai <= max_state_viol_prob;
                lb_deltai <= gammai <= max_input_viol_prob;
                % Ensure the thresholds are satisfied
                sum(deltai) <= max_state_viol_prob;
                sum(gammai) <= max_input_viol_prob;
                % Positive slack
                slack_cc_sqrt_state >= slack_cc_sqrt_state_iter/2;
                slack_cc_sqrt_input >= slack_cc_sqrt_input_iter/2;
                slack_reverse_state >= 0;
                slack_reverse_input >= 0;
                % Relaxing the norms with their slack variables
                norms(concat_input_space_A* M_matrix *...
                    sqrt_cov_concat_disturb,2,2) <= norm_input_replace_slack;
                norms(concat_safety_tube_A* (H * M_matrix + G) *...
                    sqrt_cov_concat_disturb,2,2) <= norm_state_replace_slack;
                % Enforce CC by introducing a slack variable slac_cc_sqrt_X
                %  a^T\mu - b + norm_replace * c <= slack_cc_sqrt_X^2        (a)
                %              slack_cc_sqrt_X^2 <= |m| * norm_replace*delta (b)
                % Note that (a) is a reverse convex constraint. We enforce it by
                % linearizing the RHS of (a) to obtain
                % a^T\mu - b + norm_replace * c - slack_cc_sqrt_X_iter^2 ...
                %   - 2*slack_cc_sqrt_X_iter(slack_cc_sqrt_X -
                %               slack_cc_sqrt_X_iter) <= slack_reverse_state
                % Note that (b) is a hyperbolic cone constraint which can be
                % reformulated as a second order cone constraint. 
                % (b) is true iff 
                % || [2*slack_cc_sqrt_X;      ||
                % || norm_replace - |m|*delta]||_2 <= norm_replace + |m|*delta.
                % Feasibility of these constraints implies feasibility of the
                % problem with relaxed norm constraints
                for approx_indx = 1:length(invcdf_approx_m)
                    positive_m_value = abs(invcdf_approx_m(approx_indx));
                    positive_c_value = invcdf_approx_c(approx_indx);

                    % LHS of the state CC
                    concat_safety_tube_A * mean_X - concat_safety_tube_b...
                        + norm_state_replace_slack * positive_c_value...
                        - slack_cc_sqrt_state_iter(:,approx_indx).^2 ...
                        - 2 * slack_cc_sqrt_state_iter(:,approx_indx).*...
                            (slack_cc_sqrt_state(:,approx_indx) -...
                                slack_cc_sqrt_state_iter(:,approx_indx))...
                            <= slack_reverse_state(:,approx_indx);
                    % LHS of the input CC
                    concat_input_space_A * d_vector - concat_input_space_b...
                        + norm_input_replace_slack * positive_c_value...
                        - slack_cc_sqrt_input_iter(:,approx_indx).^2 ...
                        - 2 * slack_cc_sqrt_input_iter(:,approx_indx).*...
                            (slack_cc_sqrt_input(:,approx_indx) -...
                                slack_cc_sqrt_input_iter(:,approx_indx))...
                            <= slack_reverse_input(:,approx_indx);
                    % RHS of the state CC
                    norms([2*slack_cc_sqrt_state(:,approx_indx),...
                           norm_state_replace_slack -...
                           positive_m_value * deltai],2,2) <=...
                                    norm_state_replace_slack +...
                                        positive_m_value * deltai;
                    % RHS of the input CC
                    norms([2*slack_cc_sqrt_input(:,approx_indx),...
                           norm_input_replace_slack -...
                           positive_m_value * gammai],2,2) <=...
                                    norm_input_replace_slack +...
                                        positive_m_value * gammai;
                end
        cvx_end
        
        solver_status = cvx_status;
        % Iteration status analysis
        if strcmpi(cvx_status, 'Solved') ||...
                strcmpi(cvx_status, 'Inaccurate/Solved')
            obj_curr = cvx_optval;
            dc_slack_with_tau_curr = tau_iter *...
                (sum(sum(slack_reverse_state))+sum(sum(slack_reverse_input)));    
                
            if iter_count == 0
                % Exit condition is still 1 since we want to do another
                % iteration
                if options.verbose >= 2
                    fprintf([' 0. CVX status: %s | Max iterations : %d\n',...
                             'Current norm_replace_slack val: %1.2e | ',...
                                'tau_iter: %d\n',...
                             'DC slack-total sum --- state: %1.2e | ',...
                                'input: %1.2e\n\n'],...
                            solver_status,  obj_curr, options.max_iter,...
                            tau_iter, sum(sum(slack_reverse_state)),...
                            sum(sum(slack_reverse_input)));    
                end
            else
                %The continue criteria is \leq iter_max AND 
                % NOT OF DC stopping criteria in Lipp and Boyd is met)
                continue_condition = ((iter_count <= options.iter_max) &&...
                    ~(abs((obj_prev + dc_slack_with_tau_prev)...
                           - (obj_curr + dc_slack_with_tau_curr))...
                        <= options.dc_conv_tol));
                if options.verbose >= 2
                    % Iteration status analysis
                    fprintf(['%2d. CVX status: %s | Max iterations : %d\n',...
                             'Current norm_replace_slack val: %1.2e |'...
                                ' tau_iter: %d\n',...
                             'Curr-Prev --- DC: %1.2e | Obj: %1.2e\n',...
                             'DC slack-total sum --- state: %1.2e | ',...
                                'input: %1.2e\n',...
                             'DC convergence error: %1.2e\n\n'],...
                             iter_count, solver_status, options.max_iter,...
                             obj_curr, tau_iter,... 
                             dc_slack_with_tau_prev - dc_slack_with_tau_curr,...
                             obj_prev - obj_curr,...
                             sum(sum(slack_reverse_state)),...
                             sum(sum(slack_reverse_input)),...
                             abs((obj_prev + dc_slack_with_tau_prev)...
                               - (obj_curr + dc_slack_with_tau_curr)));    
                end
            end    
            % Next iteration initialization
            slack_cc_sqrt_state_iter = slack_cc_sqrt_state;
            slack_cc_sqrt_input_iter = slack_cc_sqrt_input;            
            tau_iter = min(tau_iter * options.scaling_tau, options.tau_max);
            % Increment counter 
            iter_count = iter_count + 1;
        else
            % Converged to an infeasible solution => Quit!
            continue_condition = -1;
            if options.verbose >= 2
                disp('Converged to an infeasible solution!');
            end
        end
    end

    if sum(sum(slack_reverse_state))>= options.dc_conv_tol ||...
            sum(sum(slack_reverse_input)) >= options.dc_conv_tol ||...
            continue_condition == -1 ||...
            any(abs(norm_state_replace_slack - norms(concat_safety_tube_A *...
                (H * M_matrix + G) * sqrt_cov_concat_disturb,2,2)) >...
                    options.slack_tol) ||...
            any(abs(norm_input_replace_slack - norms(concat_input_space_A *...
                 M_matrix * sqrt_cov_concat_disturb,2,2)) > options.slack_tol)            
        dc_feas = false;                
    else    
        dc_feas = true;
    end
    sol_struct.risk_alloc_state = deltai;
    sol_struct.risk_alloc_input = gammai;
    sol_struct.stoch_reach_prob = 1-sum(deltai);
    sol_struct.input_satisfaction_prob = 1-sum(gammai);
    sol_struct.opt_input_vec = d_vector;
    sol_struct.opt_input_gain = M_matrix;
end