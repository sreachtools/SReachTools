function [lb_stoch_reach, opt_input_vec, opt_input_gain, risk_alloc_state, ...
    risk_alloc_input] = SReachPointCcA(sys, initial_state, safety_tube, options)
% Solve the problem of stochastic reachability of a target tube (a lower bound
% on the maximal reach probability and an affine controller synthesis) using
% chance-constrained optimization and difference of convex programming
% =============================================================================
%
% SReachPointCcA implements the chance-constrained underapproximation to the
% problem of stochastic reachability of a target tube to construct an affine
% controller. This technique is discussed in detail in the paper,
%
% A. Vinod and M. Oishi. Affine controller synthesis for stochastic reachability
% via difference of convex programming. In Proc. Conf. Dec. & Ctrl., 2019.
% (submitted). https://hscl.unm.edu/affinecontrollersynthesis/
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
%
% =============================================================================
%
% [lb_stoch_reach, opt_input_vec, opt_input_gain, risk_alloc_state, ...
%   risk_alloc_input] = SReachPointCcA(sys, initial_state, safety_tube, options)
% 
% Inputs:
% -------
%   sys          - System description (LtvSystem/LtiSystem object)
%   initial_state- Initial state for which the maximal reach probability must be
%                  evaluated (A numeric vector of dimension sys.state_dim)
%   safety_tube  - Collection of (potentially time-varying) safe sets that
%                  define the safe states (Tube object)
%   options      - Collection of user-specified options for 'chance-affine'
%                  (Matlab struct created using SReachPointOptions)
%
% Outputs:
% --------
%   lb_stoch_reach 
%               - Lower bound on the stochastic reachability of a target tube
%                 problem computed using chance constraints and
%                 difference-of-convex techniques. While it is expected to lie
%                 in [0,1], it is set to -1 in cases where the
%                 difference-of-convex optimization fails to converge.
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
%   risk_alloc_state 
%               - Risk allocation for the state constraints
%   risk_alloc_input
%               - Risk allocation for the input constraints
%
% See also SReachPoint.
%
% Notes:
% * We recommend using this function through SReachPoint.
% * This function requires CVX to work.
% * This function returns a **lower bound to the maximal reach probability under
%   hard input constraints**. This lower bound is obtained by a linear
%   transformation of the maximal reach probability associated with the
%   unsaturated affine controller using the user-specified likelihood threshold
%   on the hard input constraints. See Theorem 1 of the paper cited above.
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
    time_horizon = length(safety_tube) - 1;

    otherInputHandling(sys, options, time_horizon);

    % Get half space representation of the target tube and time horizon
    % skipping the first time step
    [concat_safety_tube_A, concat_safety_tube_b] = safety_tube.concat(...
        [2 time_horizon+1]);

    %% Halfspace-representation of U^N, H, G,mean_X_sans_input, cov_X_sans_input
    % GUARANTEES: Non-empty input sets (polyhedron)
    [concat_input_space_A, concat_input_space_b] = getConcatInputSpace(sys, ...
        time_horizon);
    % GUARANTEES: Compute the input concat and disturb concat transformations
    [~, H, G] = getConcatMats(sys, time_horizon);
    % GUARANTEES: well-defined initial_state and time_horizon
    sysnoi = LtvSystem('StateMatrix',sys.state_mat,'DisturbanceMatrix', ...
        sys.dist_mat,'Disturbance',sys.dist);
    X_sans_input_rv_with_init_state = SReachFwd('concat-stoch', sysnoi, ...
        initial_state, time_horizon);
    mean_X_zi_with_init_state = X_sans_input_rv_with_init_state.mean();
    mean_X_zi = mean_X_zi_with_init_state(sysnoi.state_dim + 1:end);
    
    %% Concatenation of disturbance vector
    W = sys.dist.concat(time_horizon);    
    
    %% Compute M --- the number of polytopic halfspaces to worry about
    n_lin_state = size(concat_safety_tube_A,1);
    n_lin_input = size(concat_input_space_A,1);
    
    % Compute a sparse square root of a matrix: chol produces a sqrt such that
    % sqrt' * sqrt = M. Hence, whenever, we left multiply (inside a norm), we
    % must transpose.
    [sqrt_cov_concat_disturb, p] = chol(W.cov());    
    if p > 0
        % Non-positive definite matrix can not use Cholesky's decomposition
        % Use sqrt to obtain a symmetric non-sparse square-root matrix
        sqrt_cov_concat_disturb = sqrt(W.cov());
    end
    

    %% Piecewise-affine approximation of norminvcdf
    [invcdf_approx_m, invcdf_approx_c, lb_risk] =...
        computeNormCdfInvOverApprox(0.5, options.pwa_accuracy, ...
            max(n_lin_state,n_lin_input));
        
    %% Difference of convex-based evaluation
    % Counter for the iterations        
    iter_count = 0;      

    % Initializations for DC iterative algorithm
    obj_curr = Inf;      
    % Slack initialization with gain as zero (search 'chol' for why transpose)
    norm_state_replace_slack_iter = norms(concat_safety_tube_A * G * ...
                                                sqrt_cov_concat_disturb',2,2);
    norm_input_replace_slack_iter = zeros(n_lin_input,1);
    norminvdeltai_iter = norminv(lb_risk * ones(n_lin_state,1));
    norminvgammai_iter = norminv(lb_risk * ones(n_lin_input,1));
    tau_iter = options.tau_initial;

    continue_condition = 1;    
    % DC subproblems
    while continue_condition == 1
        % Store previous iterations
        obj_prev = obj_curr;
        if options.verbose >= 2
            disp('Setting up the CVX problem');
        end
        % The iteration values are updated at the end of the problem
        cvx_begin quiet
            variable M(sys.input_dim*time_horizon,sys.dist_dim*time_horizon);
            variable d(sys.input_dim * time_horizon, 1);
            variable mean_X(sys.state_dim * time_horizon, 1);
            % State chance constraint
            variable deltai(n_lin_state, 1) nonnegative;
            variable norminvdeltai(n_lin_state, 1) nonnegative;
            variable norm_state_replace_slack(n_lin_state, 1) nonnegative;
            variable slack_reverse_state(n_lin_state, 1) nonnegative;
            % Input chance constraint
            variable gammai(n_lin_input, 1) nonnegative;
            variable norminvgammai(n_lin_input, 1) nonnegative;
            variable norm_input_replace_slack(n_lin_input, 1) nonnegative;
            variable slack_reverse_input(n_lin_input, 1) nonnegative;
            % Minimize slack variable for the norm replacements (epigraph
            % construction) and also the DC prog.-based slack constraints
            minimize (sum(deltai) + tau_iter * ...
                                        (sum(sum(slack_reverse_state)) + ...
                                            sum(sum(slack_reverse_input))));
            subject to
                % Causality constraints on M_matrix
                for time_indx = 1:time_horizon
                    M((time_indx-1)*sys.input_dim + 1:...
                        time_indx*sys.input_dim, ...
                        (time_indx-1)*sys.dist_dim+1:end) == 0; 
                end
                % slack variables
                slack_reverse_state >= 0;
                slack_reverse_input >= 0;
                norm_state_replace_slack >= 0;
                norm_input_replace_slack >= 0;
                norminvdeltai >= 0;
                norminvgammai >= 0;
                % Mean trajectory constraint (mean_X_zi already has Zx + GW)
                mean_X == mean_X_zi + H * (M * W.mean() + d);
                % Risk allocation bounds --- state
                lb_risk <= deltai <= 0.5;
                sum(deltai) <= 1 - options.max_input_viol_prob;
                % Risk allocation bounds --- input
                lb_risk <= gammai  <= options.max_input_viol_prob;
                sum(gammai) <= options.max_input_viol_prob;
                % Norms in their epigraph form (search 'chol' for why transpose)
                norms(concat_safety_tube_A* (H * M + G) * ...
                    sqrt_cov_concat_disturb',2,2)<= norm_state_replace_slack;
                norms(concat_input_space_A* M * ...
                    sqrt_cov_concat_disturb',2,2)<= norm_input_replace_slack;
                % Norminvcdf(1-x) in their epigraph form via
                % piecewise-affine approximation
                for deltai_indx = 1:n_lin_state
                    norminvdeltai(deltai_indx) >= invcdf_approx_m.* ...
                        deltai(deltai_indx) + invcdf_approx_c; 
                end
                for gammai_indx = 1:n_lin_input
                    norminvgammai(gammai_indx) >= invcdf_approx_m.* ...
                        gammai(gammai_indx) + invcdf_approx_c; 
                end

                % State CC
                concat_safety_tube_A * mean_X + ...
                  pow_p(norm_state_replace_slack + norminvdeltai,2)/2 ...
                  - concat_safety_tube_b...
                    <= norm_state_replace_slack_iter.^2/2 + ...
                            norm_state_replace_slack_iter.* ...
                                (norm_state_replace_slack - ...
                                    norm_state_replace_slack_iter) + ...
                       norminvdeltai_iter.^2/2 + norminvdeltai_iter.* ...
                                (norminvdeltai - norminvdeltai_iter) + ...
                       slack_reverse_state;
                   
                % Input CC
                concat_input_space_A * d + ...
                  pow_p(norm_input_replace_slack + norminvgammai,2)/2 ...
                  - concat_input_space_b...
                    <= norm_input_replace_slack_iter.^2/2 + ...
                            norm_input_replace_slack_iter.* ...
                                (norm_input_replace_slack - ...
                                    norm_input_replace_slack_iter) + ...
                       norminvgammai_iter.^2/2 + norminvgammai_iter.* ...
                                (norminvgammai - norminvgammai_iter) + ...
                       slack_reverse_input;
        cvx_end

        % Post solve analysis
        solver_status = cvx_status;
        sum_slack_rev_state = sum(sum(slack_reverse_state));
        sum_slack_rev_input = sum(sum(slack_reverse_input));

        if strcmpi(cvx_status, 'Solved') ||...
                strcmpi(cvx_status, 'Inaccurate/Solved')
            % Successfully solved the subproblem
            dc_slack_with_tau_curr = tau_iter * (sum_slack_rev_state + ...
                sum_slack_rev_input);    
            obj_curr = cvx_optval;

            if iter_count == 0
                if options.verbose >= 2
                    fprintf([' 0. CVX status: %s | Max iterations : <%d\n', ...
                             'Current probabilty: %1.3f | tau_iter: %d\n', ...
                             'DC slack-total sum --- state: %1.2e | ', ...
                                'input: %1.2e\n\n'], ...
                          solver_status,  options.iter_max, ...
                          1-(obj_curr - dc_slack_with_tau_curr), tau_iter, ...
                          sum_slack_rev_state,sum_slack_rev_input);    
                end
            else
                % The continue criteria is < iter_max AND 
                % NOT OF DC stopping criteria in Lipp and Boyd is met) AND
                % NOT OF slack is an acceptable replacement
                continue_condition = ((iter_count < options.iter_max) &&...
                    ~((abs(obj_prev - obj_curr) <= options.dc_conv_tol) &&...
                        max(sum_slack_rev_input, sum_slack_rev_input) <=...
                        options.slack_tol));
                if options.verbose >= 2
                    % Iteration status analysis
                    fprintf(['%2d. CVX status: %s | Max iterations : <%d\n', ...
                             'Current probabilty: %1.3f | tau_iter: %1.3e\n',...                             
                             'DC slack-total sum --- state: %1.2e | ', ...
                                'input: %1.2e | Acceptable: <%1.3e\n', ...
                             'DC convergence error: %1.2e | Acceptable:', ...
                             ' <%1.3e\n\n'], ...
                             iter_count, solver_status, options.iter_max, ...
                             1-(obj_curr - dc_slack_with_tau_curr), ...
                             tau_iter, sum_slack_rev_state, ...
                             sum_slack_rev_input, options.slack_tol, ...
                             abs(obj_prev - obj_curr), options.dc_conv_tol);
                end
            end    
            % Next iteration initialization
            norm_state_replace_slack_iter = norm_state_replace_slack;
            norminvdeltai_iter = norminvdeltai;
            norm_input_replace_slack_iter = norm_input_replace_slack;
            norminvgammai_iter = norminvgammai;
            tau_iter = min(tau_iter * options.scaling_tau, options.tau_max);
            % Increment counter 
            iter_count = iter_count + 1;
        else
            % Converged to an infeasible solution => Quit!
            continue_condition = -1;
            % Print reasons for failure
            if options.verbose >= 1
                fprintf('CVX had trouble finding solution. ', ...
                    'CVX status: %s\n', cvx_status);
                fprintf(['Slack variables of the difference-of-convex ', ...
                         'is not small enough\nDC sum-total slack --- ', ...
                         'state: %1.3e | input: %1.3e | Acceptable: ', ...
                         '<%1.1e\n'], ...
                         sum_slack_rev_state, sum_slack_rev_input, ...
                         options.dc_conv_tol);
            end
        end
    end
    
    if max(sum_slack_rev_state, sum_slack_rev_input) <= options.dc_conv_tol
        % Both the DC slack variables are below tolerance
        lb_stoch_reach = 1 - sum(deltai)/(1-options.max_input_viol_prob);
        opt_input_vec = d;
        opt_input_gain = M;
        risk_alloc_state = deltai;
        risk_alloc_input = gammai;
    else
        % Tell SReachPoint that no solution was found
        lb_stoch_reach = -1;
        opt_input_vec = nan(sys.input_dim * time_horizon,1);
        opt_input_gain = [];
        risk_alloc_state = nan(n_lin_state,1);
        risk_alloc_input = nan(n_lin_input,1);
    end
end

function otherInputHandling(sys, options, time_horizon)
    % 1. Get the correct options
    % 2. Check if the system is stochastic
    % 3. Check if the random vector is Gaussian
    % 4. Check if the disturbance matrix is identity (Theory requirement,
    %    else this formulation is not practical as it is stronger than
    %    state feedback control law)
    % NOTE: chance-affine-uni has an exact same structure
    if ~(strcmpi(options.prob_str, 'term') &&...
            strcmpi(options.method_str, 'chance-affine'))
        throwAsCaller(SrtInvalidArgsError('Invalid options provided'));
    end
    if ~isa(sys.dist,'RandomVector')
        throwAsCaller(SrtInvalidArgsError('Expected a stochastic system'));
    end
    if ~strcmpi(sys.dist.type,'Gaussian')
        throw(SrtInvalidArgsError(['SReachPointCcA requires Gaussian-',...
            'perturbed LTV/LTI system']));
    end 
    % Ensure that the disturbance matrix is identity
    err_string = ['Disturbance matrix must be an identity matrix.\nIf you ',...
        'must have the disturbance matrix, redefine the disturbance random ',...
        'vector.'];
    if sys.islti() && ~isequal(sys.dist_mat, eye(sys.dist_dim)) % Is F_k=I_k?
        throwAsCaller(SrtInvalidArgsError(err_string));
    elseif sys.isltv()
        not_equal_to_identity_flag = 0;
        for t=0:time_horizon
            if ~isequal(sys.dist_mat(t), eye(sys.dist_dim))
                not_equal_to_identity_flag = 1;
                break;
            end
        end
        if not_equal_to_identity_flag
            throwAsCaller(SrtInvalidArgsError(err_string));
        end    
    end
end
