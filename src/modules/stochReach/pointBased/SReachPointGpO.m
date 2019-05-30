function [approx_stoch_reach, opt_input_vec] = SReachPointGpO(sys, ...
    initial_state, safety_tube, options)
% Solve the problem of stochastic reachability of a target tube (a lower bound
% on the maximal reach probability and an open-loop controller synthesis) using
% Genz's algorithm and MATLAB's patternsearch (a nonlinear, derivative-free,
% constrained optimization solver)
% =============================================================================
%
% SReachPointGpO implements the Fourier transform-based underapproximation to
% the stochastic reachability of a target tube problem. The original problem was
% formulated (for the simpler problem of terminal hitting-time stochastic
% reach-avoid problem) in
%
% A. Vinod and M. Oishi, "Scalable Underapproximation for Stochastic
% Reach-Avoid Problem for High-Dimensional LTI Systems using Fourier
% Transforms," in IEEE Control Systems Letters (L-CSS), 2017.
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
%
% =============================================================================
%
%  [approx_stoch_reach, opt_input_vec] = SReachPointGpO(sys, initial_state, ...
%      safety_tube, options)
%
% Inputs:
% -------
%   sys          - System description (LtvSystem/LtiSystem object)
%   initial_state- Initial state for which the maximal reach probability must be
%                  evaluated (A numeric vector of dimension sys.state_dim)
%   safety_tube  - Collection of (potentially time-varying) safe sets that
%                  define the safe states (Tube object)
%   options      - Collection of user-specified options for 'genzps-open'
%                  (Matlab struct created using SReachPointOptions)
%
% Outputs:
% --------
%   approx_stoch_reach 
%               - Approximation of the stochastic reachability of a target tube
%                 problem | Returns -1 if patternsearch fails (exitflag < 1)
%   opt_input_vec
%               - Open-loop controller: column vector of dimension
%                 (sys.input_dim*N) x 1
%
% See also SReachPoint.
%
% Notes:
% ------
% * We recommend using this function through SReachPoint.
% * This function requires MATLAB's Global Optimization Toolbox for its
%   nonlinear solver 'patternsearch'.
% * This function requires CVX to work, since it uses SReachPointCcO for
%   initialization of MATLAB's 'patternsearch'.
% * This function uses Genz's algorithm (see in src/helperFunctions) instead of
%   MATLAB's Statistics and Machine Learning Toolbox's mvncdf to compute the
%   integral of the Gaussian over a polytope.
% * See @LtiSystem/getConcatMats for more information about the notation used.
% * This code is also used internally by SReachSetGpO (stochastic reach set
%   underapproximation via genzps-open method). There are a few adjustments done 
%   for computational purposes:
%   1. Using options.thresh, an inner min operation is used that is
%      convexity-preserving. This is an attempt to ensure that the quasi
%      Monte-Carlo simulation-driven optimization:
%      - does report a safety probability that is above prob_thresh, and
%      - does not spend too much time looking for global optimality, when a 
%        certificate of exceeding a lower bound suffices. 
%   2. After the optimization, the optimal value is reevaluated using a fresh
%      set of particles for generality.
% * A probabilistic bound on the overapproximation error between the provided
%   estimate and the true estimate may be found using Hoeffding's inequality.
% 
% ============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://sreachtools.github.io/license/
% 
%

    % INPUT HANDLING as well as necessary data creation
    optionsCc = SReachPointOptions('term', 'chance-open');
    optionsCc.desired_accuracy = options.desired_accuracy;
    [guess_approx_stoch_reach, guess_opt_input_vec, ~, extra_info] =...
        SReachPointCcO(sys, initial_state, safety_tube, optionsCc);

    otherInputHandling(sys, options);
        
    % Unpack the other necessary data from SReachPointCcO
    concat_safety_tube_A = extra_info.concat_safety_tube_A;
    concat_safety_tube_b = extra_info.concat_safety_tube_b;
    concat_input_space_A = extra_info.concat_input_space_A;
    concat_input_space_b = extra_info.concat_input_space_b;
    H = extra_info.H;
    
    % Setup for probability computation
    X_sans_input = RandomVector.gaussian(extra_info.mean_X_sans_input, ...
        extra_info.cov_X_sans_input);
    poly_tube = Polyhedron('A', concat_safety_tube_A, 'b',concat_safety_tube_b);
    
    % Time_horizon definition
    time_horizon = length(safety_tube) - 1;
    
    if guess_approx_stoch_reach < 0
        % Chance constrained approach failed to find an optimal open-loop
        % controller => Try just getting the mean to be as safe as possible
        % 
        % minimize ||s||_1
        % subject to 
        %                               mean_X = mean_X_sans_input + H U
        %  concatenated_safety_tube.A * mean_X <= concatenated_safety_tube.b + s   
        %                                                (Softened reach const.)
        %       concatenated_input_space.A * U <= concatenated_input_space.b
        %
        % Here, mean_X_sans_input is concatenated_A_matrix * x_0 + G_matrix *
        % concatenated_mean_vector_for_the_disturbance (see 
        % @LtiSystem/getConcatMats for more information on the notation)
        %
        % OVERWRITE the guess_opt_input_vec solution from SReachPointCcO
        cvx_begin quiet
            variable U(sys.input_dim * time_horizon,1);
            variable mean_X(sys.state_dim * time_horizon,1);
            variable slack_vars(length(concat_safety_tube_b)) nonnegative;
            minimize sum(slack_vars)
            subject to
                mean_X == X_sans_input.mean() + H * U;
                concat_safety_tube_A * mean_X <= concat_safety_tube_b...
                                                    + slack_vars;
                concat_input_space_A * U <= concat_input_space_b;
        cvx_end
        if strcmpi(cvx_status,'Solved')
            % We have an initial solution to work with
            guess_approx_stoch_reach = 0;
        else
            % This should never be happening for a non-empty since it is a slack
            % variable relaxed problem
        end
    end

    if guess_approx_stoch_reach < 0
        %% No initial solution to work with
        approx_stoch_reach = -1;
        opt_input_vec = nan(sys.input_dim * time_horizon, 1);
    else
        % Minimum saturated ReachProb(U)=\int_{safety_tube}\psi_{\bX}(Z;x_0,U)dZ
        % Specifically, we saturate the reach probability from below by 
        % options.desired_accuracy => max(options.desired_accuracy,ReachProb(U))
        % in order to be able to log | RandomVector/getProbPolyhedron returns an
        % underapproximation of the reach probability
        MinSatReachProbGivenInputVec = @(input_vec) max( ...
            options.desired_accuracy, ...
            getProbPolyhedron(X_sans_input + H * input_vec, poly_tube, ...
                options.desired_accuracy));
        % Construct the reach cost function
        % Default value for options.thresh is 1
        if 1 - options.thresh > options.desired_accuracy
            if options.thresh < options.desired_accuracy
                % Complain
                throw(SrtInvalidArgsError(['Threshold can not be smaller ', ...
                    'than the desired accuracy']));
            else
                % For options.thresh < 1, we do not care about feasibility beyond
                % options.thresh. Therefore, we maximize
                %           min(ReachProb(U), options.thresh)
                % instead of ReachProb(U). Therefore, we wish to minimize
                %          -log(min(ReachProb(U), options.thresh)).
                negLogReachProbGivenInputVec= @(input_vec) -log( min( ...
                    options.thresh, MinSatReachProbGivenInputVec(input_vec)));
            end
        else
            % We wish to minimize -log(ReachProb(U))
            negLogReachProbGivenInputVec= @(input_vec) ...
                -log(MinSatReachProbGivenInputVec(input_vec));
        end

        % Compute the optimal admissible input_vec that minimizes
        % negLogReachProbGivenInputVec(input_vec) given x_0 | We use
        % patternsearch because negLogReachProbGivenInputVec can be noisy
        % due to its reliance on Monte Carlo simulation
        %
        % minimize -log(ReachProb(U))
        % subject to
        %        concat_input_space_A * U <= concat_input_space_b
        [opt_input_vec,~, exitflag] = patternsearch(...
            negLogReachProbGivenInputVec, guess_opt_input_vec, ...
            concat_input_space_A, concat_input_space_b, [],[],[],[],[], ...
            options.PSoptions);
        
        % Compute the lower bound and the optimal open_loop_control_policy
        % using a fresh set of samples
        if exitflag >= 1
            approx_stoch_reach = exp( ...
                -negLogReachProbGivenInputVec(opt_input_vec));
        else
            approx_stoch_reach = -1;
            opt_input_vec = nan(sys.input_dim * time_horizon, 1);
        end
    end
end

function otherInputHandling(sys, options)
    % 1. Get the correct options
    % 2. Check if the system is stochastic
    % 3. Check if the random vector is Gaussian    
    if ~(strcmpi(options.prob_str, 'term') &&...
            strcmpi(options.method_str, 'genzps-open'))
        throwAsCaller(SrtInvalidArgsError('Invalid options provided'));
    end
    if ~isa(sys.dist,'RandomVector')
        throwAsCaller(SrtInvalidArgsError('Expected a stochastic system'));
    end
    if ~strcmpi(sys.dist.type,'Gaussian')
        throw(SrtInvalidArgsError(['SReachPointGpO requires Gaussian-',...
            'perturbed LTV/LTI system']));
    end
end
