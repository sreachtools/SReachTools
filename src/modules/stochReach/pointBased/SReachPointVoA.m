function [approx_stoch_reach, opt_input_vec, opt_input_gain] = SReachPointVoA(...
    sys, initial_state, safety_tube, options)
% Solve the problem of stochastic reachability of a target tube (a lower bound
% on the maximal reach probability and an affine controller synthesis) using
% under-sampled particle control and coordinate descent algorithm
% =============================================================================
%
% SReachPointVoA implements the under-sampled particle control to the problem of
% stochastic reachability of a target tube to construct an affine controller.
% This technique is discussed in detail in the paper,
%
% TODO: Add paper
% TODO: Add summary
%
% =============================================================================
%
%   [lb_stoch_reach, opt_input_vec, opt_input_gain, ...
%    risk_alloc_state, risk_alloc_input] = SReachPointVoA(sys,...
%       initial_state, safety_tube, options)
% 
% Inputs:
% -------
%   sys          - System description (LtvSystem/LtiSystem object)
%   initial_state- Initial state for which the maximal reach probability must be
%                  evaluated (A numeric vector of dimension sys.state_dim)
%   safety_tube  - Collection of (potentially time-varying) safe sets that
%                  define the safe states (Tube object)
%   options      - Collection of user-specified options for 'voronoi-affine'
%                  (Matlab struct created using SReachPointOptions)
%
% Outputs:
% --------
%   lb_stoch_reach 
%               - Lower bound on the stochastic reachability of a target tube
%                 problem computed using chance constraints and
%                 difference-of-convex techniques
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
    opt_input_gain = [];

    % Requires Gurobi since we are solving a MILP
    [default_solver, solvers_cvx] = cvx_solver;
    if ~(contains(default_solver,'Gurobi') ||...
        any(contains(solvers_cvx,'Gurobi')))
        warning('SReachTools:setup',['SReachPointVoA returns a trivial ', ...
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
        n_lin_input = size(concat_input_space_A,1);

        % Step 1: Compute the number of particles needed to meet the
        % specified tolerance
        n_particles = ceil(-log(options.failure_risk)/...
            (2*options.max_overapprox_err));        
        n_kmeans = max(ceil(n_particles * options.undersampling_fraction),...
            options.min_samples);
        if n_kmeans > 100
            warning('SReachTools:runtime',sprintf(['Particle control with ',...
                'more than 100 samples may cause computational problems.\n',...
                'Going to analyze %d samples.'], n_kmeans));
        elseif n_kmeans == options.min_samples && options.verbose >= 1
            disp(['Undersampling fraction is too severe. Use the minimum',...
                 ' number of particles prescribed.']);
        end
        
        if options.verbose >= 1
            fprintf(['Required number of particles: %1.4e | Samples ',...
                'used: %3d\n'], n_particles, n_kmeans);
            fprintf('Creating Gaussian random variable realizations....');
        end        
        % Compute the stochasticity of the concatenated disturbance random vec
        W = concat(sys.dist, time_horizon);        
        % Create realizations of W arranged columnwise
        W_realizations = mvnrnd(W.parameters.mean', W.parameters.covariance,...
            n_particles)';
        if options.verbose >= 1
            fprintf('Done\n');
        end

        % Step 2: Compute the Voronoi centers with prescribed number of
        % bins --- Use k-means clustering to undersample (H*M + G)*W
        % Transposed input since kmeans expects each data point row-wise
        if options.verbose >= 1
            fprintf('Using k-means for undersampling....');
        end        
        [idx, W_centroids_output] = kmeans(W_realizations', n_kmeans,...
            'MaxIter',1000);        
        W_centroids = W_centroids_output';
        % Count the number of disturbance realizations in each of the partition
        voronoi_count = zeros(n_kmeans,1);
        for idx_indx = 1:n_kmeans
            voronoi_count(idx_indx) = nnz(idx==idx_indx);
        end
        
        % Step 3: Solve the undersampled MILP that computes buffers online
        if options.verbose >= 1
            fprintf('Done\n');
            fprintf('Setting up CVX problem....');
        end
        cvx_begin
            if options.verbose >= 2
                cvx_quiet false
            else
                cvx_quiet true
            end
            cvx_solver Gurobi
            variable M(sys.input_dim * time_horizon, sys.dist_dim*time_horizon);
            variable D(sys.input_dim * time_horizon,1);
            variable X_realization(sys.state_dim * time_horizon, n_kmeans);
            variable U_realization(sys.input_dim * time_horizon, n_kmeans);
            variable zx(1,n_kmeans) binary;
            variable zu(1,n_kmeans) binary;
            variable zxu(1,n_kmeans) binary;
            variable buffers_x(n_lin_state, n_kmeans);
            variable buffers_u(n_lin_input, n_kmeans);
            maximize ((zxu*voronoi_count)/n_particles)
            subject to
                X_realization == repmat(Z * initial_state + H * D, ...
                    1, n_kmeans) + (H * M + G) * W_centroids;
                U_realization == M * W_centroids + repmat(D, 1, n_kmeans);
                % Chance constraints for the state: Definition
                concat_safety_tube_A * X_realization + buffers_x <= repmat( ...
                    concat_safety_tube_b, 1, n_kmeans) + ...
                    options.bigM * repmat(1-zx, n_lin_state, 1);
                % Chance constraints for the input: Definition
                concat_input_space_A * U_realization + buffers_u <= repmat( ...
                    concat_input_space_b, 1, n_kmeans) + ...
                    options.bigM * repmat(1-zu, n_lin_input, 1);
                % Chance constraints for the input: Constraint imposition
                (zu*voronoi_count)/n_particles >= 1 - options.max_input_viol_prob;
                % zxu = zx && zu
                % http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.42.7380&rep=rep1&type=pdf
                zxu <= zx;
                zxu <= zu;
                zx + zu <= 1 + zxu;
                for idx_indx = 1:n_kmeans
                    % buffers_x: Definition
                    buffers_x(:, idx_indx) >= max(concat_safety_tube_A *...
                        (H * M + G) * (W_realizations(:, idx==idx_indx) -...
                            W_centroids(:, idx_indx)), [], 2);
                    % buffers_u: Definition                
                    buffers_u(:, idx_indx) >= max(concat_input_space_A *...
                        M * (W_realizations(:, idx==idx_indx) - ...
                            W_centroids(:, idx_indx)), [], 2);
                end
            if options.verbose >= 1
                fprintf('Done\nParsing and solving the MILP....');
            end
        cvx_end       
        if options.verbose >= 1
            fprintf('Done\n');
        end
        %% Overwrite the solutions
        if strcmpi(cvx_status, 'Solved')
            approx_voronoi_stoch_reach = cvx_optval;
            opt_input_vec = D; 
            opt_input_gain = M;
            % Step 4b: Improve upon the estimate by a refined counting
            % Solve the mixed-integer linear program
            opt_U_realizations = M * W_realizations + repmat(D, 1, n_particles);
            opt_X_realizations = repmat(Z*initial_state,1,n_particles) + ...
                H * opt_U_realizations + G * W_realizations;
            % All by default does columnwise (A particle succeeds only if
            % it clears all hyperplanes --- has a column of ones)
            zx_orig = all(concat_safety_tube_A * opt_X_realizations <=...
                concat_safety_tube_b);
            zu_orig = all(concat_input_space_A * opt_U_realizations <=...
                concat_input_space_b);
            approx_stoch_reach = sum(zx_orig.*zu_orig)/n_particles;            
            if options.verbose >= 1
                fprintf(...
                    ['Undersampled probability (with %d particles): %1.3f\n',...
                    'Underapproximation to the original MILP (with %d ',...
                    'particles): %1.3f\n'],...
                    n_kmeans, approx_voronoi_stoch_reach,...
                    n_particles, approx_stoch_reach);
            end
        end
    end
end

function otherInputHandling(options, sys)
    if ~(strcmpi(options.prob_str, 'term') &&...
            strcmpi(options.method_str, 'voronoi-affine'))
        throwAsCaller(SrtInvalidArgsError('Invalid options provided'));
    end
    if ~strcmpi(sys.dist.type,'Gaussian')
        throwAsCaller(SrtInvalidArgsError(['Expected a Gaussian-perturbed', ...
            'linear system']));
    end
end

% Implemented the buffer assignment via max(M*(W^{(j)}_i - W^{(j)}_c) <=
% buffer which in turn is equivalent to a collection of linear constraints.
% Currently implements the chance constraint approach
% Cost is zx & zu since we need P{X safe | acceptable input}.
% 1. Fix comments
% 2. Fix input handling
% 3. Fix docstrings here and in SReachPointOptions