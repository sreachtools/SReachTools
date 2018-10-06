function [stoch_reach_prob, opt_input_vec, opt_input_gain,...
    risk_alloc_state, risk_alloc_input] = SReachPointCcA(sys, initial_state,...
        safety_tube, options)
% Solve the stochastic reach-avoid problem (lower bound on the probability and
% an affine controller synthesis) using chance-constrained convex optimization
% =============================================================================
%
% SReachPointCcA implements the chance-constrained convex underapproximation to 
% the problem of stochastic reachability of a target tube
%
% A. Vinod and M. Oishi, HSCC, 2019 (TODO)
%
% This function uses difference-of-convex algorithm (also known as the 
% convex-concave procedure) to compute a local optima for the risk
% allocation and an associated affine controller.
%
% =============================================================================
%
%   [stoch_reach_prob, opt_input_vec, opt_input_gain,...
%       risk_alloc_state, risk_alloc_input] = SReachPointCcA(sys,...
%        initial_state, safety_tube, options)
% 
% Inputs:
% -------
%   sys          - System description (LtvSystem/LtiSystem object)
%   initial_state- Initial state for which the maximal reach probability must be
%                  evaluated (A numeric vector of dimension sys.state_dim)
%   safety_tube  - Collection of (potentially time-varying) safe sets that
%                  define the safe states (TargetTube object)
%   options      - Collection of user-specified options for 'chance-affine'
%                  (Matlab struct created using SReachPointOptions)
%
% Outputs:
% --------
%   stoch_reach_prob 
%               - Lower bound on the stochastic reachability of a target 
%                 tube problem computed using convex chance
%                      constraints and difference-of-convex techniques
%   opt_input_vec, 
%     opt_input_gain
%               - Controller U=MW+d for a concatenated input vector 
%                   U = [u_0; u_1; ...; u_{N-1}] and concatenated disturbance
%                   vector W=[w_0; w_1; ...; w_{N-1}]. 
%                   - opt_input_gain: Affine controller gain matrix of dimension
%                       (sys.input_dim*N) x (sys.dist_dim*N)
%                   - opt_input_vec: Open-loop controller: column vector dimension
%                       (sys.input_dim*N) x 1
%   risk_alloc_state 
%               - Risk allocation for the state constraints
%   risk_alloc_input
%               - Risk allocation for the input constraints
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
    inpar.addRequired('initial_state', @(x) validateattributes(x,...
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
    time_horizon = length(safety_tube) - 1;

    % Get half space representation of the target tube and time horizon
    % skipping the first time step
    [concat_safety_tube_A, concat_safety_tube_b] = safety_tube.concat(...
        [2 time_horizon+1]);

    %% Halfspace-representation of U^N, H, G,mean_X_sans_input, cov_X_sans_input
    % GUARANTEES: Non-empty input sets (polyhedron)
    [concat_input_space_A, concat_input_space_b] = getConcatInputSpace(sys,...
        time_horizon);
    % GUARANTEES: Compute the input concat and disturb concat transformations
    [~, H, G] = getConcatMats(sys, time_horizon);
    % GUARANTEES: Gaussian-perturbed LTI system (sys) and well-defined
    % initial_state and time_horizon
    sysnoi = LtvSystem('StateMatrix',sys.state_mat,'DisturbanceMatrix',...
        sys.dist_mat,'Disturbance',sys.dist);
    [mean_X_sans_input, ~] = SReachFwd('concat-stoch', sysnoi, initial_state,...
        time_horizon);
    
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
            computeNormCdfInvOverApprox(...
                max(max_state_viol_prob, options.max_input_viol_prob),...
                options.pwa_accuracy,...
                max(n_lin_state,n_lin_input));
        
        %% Difference of convex-based evaluation
        % Counter for the iterations        
        iter_count = 0;      

        % Initializations for DC iterative algorithm
        obj_curr = Inf;      
        dc_slack_with_tau_curr = Inf;

        % Initialization of the difference of convex program
        n_pwa = length(invcdf_approx_m);
        slack_cc_sqrt_state_iter = zeros(n_lin_state, n_pwa);
        for approx_indx = 1:length(invcdf_approx_m)
            positive_c_value = invcdf_approx_c(approx_indx);
            % a^T * \mu - b + c|| || <= s := c|| ||
            slack_cc_sqrt_state_iter(:,approx_indx) = ...
                positive_c_value .* norms(concat_safety_tube_A * G *...
                        sqrt_cov_concat_disturb,2,2);
        end

        slack_cc_sqrt_input_iter = zeros(n_lin_input, n_pwa);            
        tau_iter = options.tau_initial;

        continue_condition = 1;    
        % DC subproblems
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
                % construction) and also the DC prog.-based slack constraints
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
                    lb_deltai <= gammai <= options.max_input_viol_prob;
                    % Ensure the thresholds are satisfied
                    sum(deltai) <= max_state_viol_prob;
                    sum(gammai) <= options.max_input_viol_prob;
                    % Positive slack
                    slack_cc_sqrt_state >= slack_cc_sqrt_state_iter/2;
                    slack_cc_sqrt_input >= slack_cc_sqrt_input_iter/2;
                    slack_reverse_state >= 0;
                    slack_reverse_input >= 0;
                    % Relaxing the norms with their slack variables
                    norms(concat_input_space_A* M_matrix *...
                        sqrt_cov_concat_disturb,2,2)<= norm_input_replace_slack;
                    norms(concat_safety_tube_A* (H * M_matrix + G) *...
                        sqrt_cov_concat_disturb,2,2)<= norm_state_replace_slack;
                    % Enforce CC by introducing a slack variable slac_cc_sqrt_X
                    %  a^T\mu-b +norm_replace*c <= slack_cc_sqrt_X^2        (a)
                    %         slack_cc_sqrt_X^2 <= |m| * norm_replace*delta (b)
                    % Note that (a) is a reverse convex constraint. We enforce 
                    % it by linearizing the RHS of (a) to obtain
                    % a^T\mu - b + norm_replace * c - slack_cc_sqrt_X_iter^2 ...
                    %   - 2*slack_cc_sqrt_X_iter(slack_cc_sqrt_X -
                    %               slack_cc_sqrt_X_iter) <= slack_reverse_state
                    % Note that (b) is a hyperbolic cone constraint which can be
                    % reformulated as a second order cone constraint. 
                    % (b) is true iff 
                    % || [2*slack_cc_sqrt_X;      ||
                    % || norm_replace - |m|*delta]||_2<=norm_replace + |m|*delta
                    % Feasibility of these constraints implies feasibility of 
                    % the problem with relaxed norm constraints
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
                        concat_input_space_A * d_vector -concat_input_space_b...
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
            
            % Post solve analysis
            solver_status = cvx_status;
            sum_slack_rev_state = sum(sum(slack_reverse_state));
            sum_slack_rev_input = sum(sum(slack_reverse_input));
            
            if strcmpi(cvx_status, 'Solved') ||...
                    strcmpi(cvx_status, 'Inaccurate/Solved')
                
                % Successfully solved the subproblem
                dc_slack_with_tau_curr = tau_iter * (sum_slack_rev_state +...
                    sum_slack_rev_input);    
                obj_curr = cvx_optval;
                
                if iter_count == 0
                    % Exit condition is still 1 since we want to do another
                    % iteration
                    if options.verbose >= 2
                        fprintf([' 0. CVX status: %s | Max iterations : <%d',...
                                 '\nCurrent sum_norm_replace_slack: %1.2e |',...
                                    ' tau_iter: %d\n',...
                                 'DC slack-total sum --- state: %1.2e | ',...
                                    'input: %1.2e\n\n'],...
                              solver_status,  options.iter_max,...
                              obj_curr - dc_slack_with_tau_curr, tau_iter,...
                              sum_slack_rev_state,sum_slack_rev_input);    
                    end
                else
                    norm_state_gain = norms(concat_safety_tube_A *...
                              (H * M_matrix + G) * sqrt_cov_concat_disturb,2,2);
                    norm_input_gain = norms(concat_input_space_A *...
                               M_matrix * sqrt_cov_concat_disturb,2,2);
                    slack_equal_norm = (max(norm_state_replace_slack -...
                                    norm_state_gain) < options.slack_tol) &&...
                                       (max(norm_input_replace_slack -...
                                    norm_input_gain) < options.slack_tol);
                    %The continue criteria is < iter_max AND 
                    % NOT OF DC stopping criteria in Lipp and Boyd is met) AND
                    % NOT OF slack is an acceptable replacement
                    continue_condition = ((iter_count < options.iter_max) &&...
                        ~((abs(obj_prev - obj_curr) <= options.dc_conv_tol)...
                            && slack_equal_norm));
                    if options.verbose >= 2
                        % Iteration status analysis
                        fprintf(['%2d. CVX status: %s | Max iterations : ',...
                                 '<%d\nCurrent sum_norm_replace_slack: %1.2e'...
                                    ' | tau_iter: %d\n',...
                                 'DC slack-total sum --- state: %1.2e | ',...
                                    'input: %1.2e\n',...
                                 'DC convergence error: %1.2e | Acceptable:',...
                                 ' <%1.3e\n\n'],...
                                 iter_count, solver_status, options.iter_max,...
                                 obj_curr - dc_slack_with_tau_curr, tau_iter,... 
                                 sum_slack_rev_state, sum_slack_rev_input,...
                                 abs(obj_prev - obj_curr), options.dc_conv_tol);    
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
                    fprintf(['CVX had trouble solving this subproblem | ',...
                        'CVX status: %s\n'],...
                        cvx_status);
                end
            end
        end
        % Post DC analysis
        if sum_slack_rev_state >= options.dc_conv_tol ||...
                sum_slack_rev_input >= options.dc_conv_tol ||...
                continue_condition == -1 || ~slack_equal_norm                        
            % Print reasons for failure
            if options.verbose >= 1
                if slack_equal_norm
                    fprintf(['Slack variables of the difference-of-convex ',...
                             'is not small enough\nDC sum-total slack --- ',...
                             'state: %1.3e | input: %1.3e | Acceptable: ',...
                             '<%1.1e\n'],...
                             sum_slack_rev_state, sum_slack_rev_input,...
                             options.dc_conv_tol);
                else
                    fprintf(['Increasing the maximum number of iterations ',...
                             'might help! Slack variables is not equal to ',...
                             'norms.\nMax error in slack --- state: %1.3e |',...
                             ' input: %1.3e | Acceptable: <%1.1e\nMin error',...
                             ' in slack --- state: %1.3e | input: %1.3e ',...
                             '(both must be +ve)\n'],...
                             max(norm_state_replace_slack - norm_state_gain),...
                             max(norm_input_replace_slack - norm_input_gain),...
                             options.slack_tol,...
                             min(norm_state_replace_slack - norm_state_gain),...
                             min(norm_input_replace_slack - norm_input_gain));                
                end
            end
            % If infeasible for 0.5 itself, then no hope!
            if abs(max_state_viol_prob - 0.5) < eps
                % Exit from the while loop by setting the gap to be zero
                ub_max_state_viol_prob = lb_max_state_viol_prob;
                % Tell SReachPoint that no solution was found
                stoch_reach_prob = -1;
                opt_input_vec = nan(sys.input_dim * time_horizon,1);
                opt_input_gain = [];
                risk_alloc_state = nan(n_lin_state,1);
                risk_alloc_input = nan(n_lin_input,1);
                if options.verbose >= 1
                    disp('Infeasible for Delta_x = 0.5! No solution expected.');
                end            
            else
                % Increase the lower bound
                lb_max_state_viol_prob = max_state_viol_prob;
                if options.verbose >= 1
                    disp(['Infeasible! Increasing Delta_x ',...
                          '(relaxing state violation probability threshold)']);
                end
            end
        else    
            % Decrease the upper bound
            ub_max_state_viol_prob = max_state_viol_prob;
            % Update the solution
            stoch_reach_prob = 1 - sum(deltai);
            opt_input_vec = d_vector;
            opt_input_gain = M_matrix;
            risk_alloc_state = deltai;
            risk_alloc_input = gammai;
            if options.verbose >= 1
                disp(['Feasible solution found! Decreasing Delta_x ',...
                      '(tightening state violation probability threshold)']);
            end
        end        
    end    
end

% %% Difference of convex subproblem implemented in CVX
%         [sol_struct_temp, dc_feas] = differenceOfConvexApproach(...
%             max_state_viol_prob,...
%             options.max_input_viol_prob,... 
%             sys,... 
%             time_horizon,...
%             n_lin_state,...
%             n_lin_input,... 
%             mean_X_sans_input,...
%             H,...
%             G,...
%             concat_safety_tube_A,...
%             concat_safety_tube_b,...
%             concat_input_space_A,...
%             concat_input_space_b,...
%             sqrt_cov_concat_disturb,...
%             invcdf_approx_m,...
%             invcdf_approx_c,...
%             lb_deltai,...
%             options);
function otherInputHandling(options)
    if ~(strcmpi(options.prob_str, 'term') &&...
            strcmpi(options.method_str, 'chance-affine'))
        throwAsCaller(SrtInvalidArgsError('Invalid options provided'));
    end
end
