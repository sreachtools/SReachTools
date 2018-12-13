function [approx_stoch_reach, opt_input_vec, opt_input_gain, kmeans_info] =...
    SReachPointVoA(sys, initial_state, safety_tube, options)
% Solve the problem of stochastic reachability of a target tube (a lower bound
% on the maximal reach probability and an affine controller synthesis) using
% under-sampled particle control and coordinate descent algorithm
% =============================================================================
%
% SReachPointVoA implements the under-sampled particle control to the problem of
% stochastic reachability of a target tube to construct an affine controller.
%
%    High-level desc.   : Sample particles based on the additive noise and solve
%                         a mixed-integer linear program to make the maximum
%                         number of particles satisfy the reachability objective
%                         In addition, we use Voronoi partition to
%                         drastically improve the tractability while
%                         preserving the underapproximation quality
%    Approximation      : Overapproximation bounded above (in probability) by a
%                         user-specified tolerance
%    Controller type    : A history-dependent affine controller that satisfies
%                         softened input constraints (controller satisfies the
%                         hard input bounds upto a user-specified probabilistic
%                         threshold)
%    Optimality         : Suboptimal (w.r.t particles drawn) affine disturbance
%                         feedback controller 
%    Dependency (EXT)   : CVX, Gurobi
%    SReachTool function: SReachPointVoA
%    Paper              : TODO
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
%   kmeans_info - A MATLAB struct containing the information about partitioning
%                 of W space. The struct contains the following info:
%                  n_particles    - Number of particles based off Hoeffding's
%                                   inequality
%                  n_kmeans       - Number of bins for kmeans clustering
%                  W_centroids    - Centroids obtained from kmeans clustering
%                  W_realizations - Realizations for the random vector W
%
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
% * Due to numerical issues, we add a small positive perturbation to the
%   b term, whenever determining containment --- Ax<=b
% * See @LtiSystem/getConcatMats for more information about the notation used.
% 
% ============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 
%
    
    % Due to numerical issues, we add a small positive perturbation to the
    % b term whenever determining containment --- Ax<=b
    my_eps = 1e-10;
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

    % Ensure that system is stochastic
    if ~isa(sys.dist,'RandomVector')
        throwAsCaller(SrtInvalidArgsError('Expected a stochastic system'));
    end
    
    % Target tubes has polyhedra T_0, T_1, ..., T_{time_horizon}
    time_horizon = length(safety_tube)-1;

    approx_stoch_reach = -1;
    opt_input_vec = nan(sys.input_dim * time_horizon,1);
    opt_input_gain = [];
    kmeans_info = [];

    % Requires Gurobi since we are solving a MILP
    [default_solver, solvers_cvx] = cvx_solver;
    if ~(contains(default_solver,'Gurobi') ||...
        any(contains(solvers_cvx,'Gurobi')))
        warning('SReachTools:setup',['SReachPointVoA returns a trivial ', ...
            'result since Gurobi (a MILP solver) was not setup.']);
    else
        % Ensure options is good
        otherInputHandling(options);
        
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
        n_particles = ceil(...
           git (log(2)-log(options.failure_risk))/(2*options.max_overapprox_err^2));        
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

        % Step 2: Obtain the random vector realizations
        if options.verbose >= 1
            fprintf(['Required number of particles: %1.4e | Samples ',...
                'used: %3d\n'], n_particles, n_kmeans);
            fprintf('Creating random variable realizations....');
        end        
        % Compute the stochasticity of the concatenated disturbance random vec
        W = concat(sys.dist, time_horizon);        
        % Create realizations of W arranged columnwise
        W_realizations = W.getRealizations(n_particles);
        if options.verbose >= 1
            fprintf('Done\n');
        end

        % Step 3: Compute the Voronoi centers with prescribed number of bins
        % --- Use k-means clustering to undersample W | Count the particles
        % associated with each center | Transposed input since kmeans expects
        % each data point row-wise
        
        if n_kmeans < n_particles
            % No. of centroids required is smaller than the actual number
            % of seeds
            if options.verbose >= 1
                fprintf('Using k-means for undersampling....');
            end 
            [idx, W_centroids_output] = kmeans(W_realizations', n_kmeans,...
                'MaxIter',1000);        
            W_centroids = W_centroids_output';
            voronoi_count = zeros(n_kmeans,1);
            for idx_indx = 1:n_kmeans
                relv_idx = (idx==idx_indx);
                % Count the particles associated with the partition
                voronoi_count(idx_indx) = nnz(relv_idx);
            end
        else
            warning('SReachTools:runtime',['Skipping Voronoi',...
                ' partitioning, since no. of centroids is equal to no. of ',...
                'particles | However, the code used will be the same.']);
            voronoi_count = ones(n_particles,1);
            W_centroids = W_realizations;     
            idx = 1:n_particles;
        end
        
        
        % Step 4a: Solve the undersampled MILP that computes buffers online
        if options.verbose >= 1
            fprintf('Done\n');
            fprintf('Setting up CVX problem....\n');
            % Suppress CVX warnings, if any
            evalc('cvx_begin;cvx_end');
            fprintf(['Fraction of buffer constraints enforced (out of 1):',...
                '0.000\n']);            
        end
        % Normalize the hyperplanes so that Ax <= b + M is good. Here, b<=1e-3
        % or 1
        [concat_safety_tube_A, concat_safety_tube_b] =...
            normalizeForParticleControl(concat_safety_tube_A,...
                concat_safety_tube_b);
        [concat_input_space_A, concat_input_space_b] =...
            normalizeForParticleControl(concat_input_space_A,...
                concat_input_space_b);
            
        % No. of non-zero blocks
        n_blocks = (time_horizon - 1) * time_horizon/2;
        blocks_indx_vec = cumsum(1:time_horizon-1);
        % CVX problem setup
        cvx_clear
        cvx_begin
            if options.verbose >= 2
                cvx_quiet false
            else
                cvx_quiet true
            end
            cvx_solver Gurobi
            expression M(sys.input_dim*time_horizon, sys.dist_dim*time_horizon);
            variable M_vars(sys.input_dim * sys.dist_dim, n_blocks);
            variable D(sys.input_dim * time_horizon,1);
            variable X_centroids(sys.state_dim * time_horizon, n_kmeans);
            variable U_centroids(sys.input_dim * time_horizon, n_kmeans);
            variable bin_x(1,n_kmeans) binary;
            variable bin_u(1,n_kmeans) binary;
            maximize ((bin_x * voronoi_count)/n_particles)
            subject to
                % Causality constraints on M is automatically enforced by
                % expression declaration (sets all other terms to zero)
                for time_indx = 2:time_horizon
                    if time_indx == 2
                        blocks_start_indx = 1;
                    else
                        blocks_start_indx = blocks_indx_vec(time_indx - 2)+1;
                    end
                    blocks_end_indx = blocks_indx_vec(time_indx - 1);
                    M((time_indx-1)*sys.input_dim + 1:...
                        time_indx*sys.input_dim, ...
                        1:(time_indx-1)*sys.dist_dim) =...
                        reshape(M_vars(:,blocks_start_indx:blocks_end_indx),...
                            sys.input_dim, []); 
                end
                % X = Z x_0 + H D + (H M + G) W
                X_centroids == repmat(Z*initial_state + H*D, 1, n_kmeans)+ ...
                    (H * M + G) * W_centroids;
                % U = D + M W
                U_centroids == M * W_centroids + repmat(D, 1, n_kmeans);
                % Chance constraints for the input: Constraint imposition
                (bin_u * voronoi_count)/n_particles >= ...
                    1 - options.max_input_viol_prob +...
                    options.max_overapprox_err;
                % Chance constraints: Definition
                for idx_indx = 1:n_kmeans
                    if options.verbose >= 1
                        fprintf('\b\b\b\b\b\b%1.3f\n', idx_indx/n_kmeans);
                    end
                    % Displacement of actual realizations from the centroids
                    W_disp = W_realizations(:, idx==idx_indx) -...
                        W_centroids(:,idx_indx);
                    % Chance constraints for the state centroids: Definition
                    % concat_safety_tube_A * X_centroid +
                    %   concat_safety_tube_A * (HM + G) * W_disp <=
                    %       concat_safety_tube_b + bigM ( 1-bin_x)
                    concat_safety_tube_A * X_centroids(:, idx_indx) -...
                        concat_safety_tube_b-options.bigM*(1-bin_x(idx_indx))+...
                        max(concat_safety_tube_A *(H * M + G) * W_disp,[],2)<=0;
                    % Chance constraints for the input centroids: Definition
                    % concat_input_space_A * U_centroid +
                    %   concat_input_space_A * M * W_disp <=
                    %       concat_input_space_b + bigM ( 1-bin_u)
                    concat_input_space_A * U_centroids(:, idx_indx) -...
                        concat_input_space_b-options.bigM*(1-bin_u(idx_indx))+...
                        max(concat_input_space_A * M * W_disp,[],2)<=0;
                end
            if options.verbose >= 1
                disp('Setup of CVX problem complete');
                fprintf('Parsing and solving the MILP....');
            end
        cvx_end       
        if options.verbose >= 1
            fprintf('Done\n');
        end
        %% Overwrite the solutions
        switch cvx_status
            case {'Solved','Inaccurate/Solved'}
                approx_voronoi_stoch_reach = cvx_optval;
                opt_input_vec = D; 
                opt_input_gain = M;
                % Step 4b: Improve upon the estimate by a refined counting
                % Solve the mixed-integer linear program
                opt_U_realizations = M * W_realizations +...
                    repmat(D, 1, n_particles);
                opt_X_realizations = repmat(Z*initial_state,1,n_particles) + ...
                    H * opt_U_realizations + G * W_realizations;
                
                % All by default does columnwise (A particle succeeds only if
                % it clears all hyperplanes --- has a column of ones) | Due to 
                % numerical issues, we add a small positive perturbation to the
                % b term, whenever determining containment --- Ax<=b
                bin_x_orig = all(concat_safety_tube_A * opt_X_realizations <=...
                    concat_safety_tube_b + my_eps);
                approx_stoch_reach_noc =sum(bin_x_orig)/n_particles;            
                
                % Step 5: Account for the saturation via HSCC 2019 result
                approx_stoch_reach = (approx_stoch_reach_noc -...
                   options.max_input_viol_prob)/(1-options.max_input_viol_prob);
                if options.verbose >= 1
                    bin_u_orig = all(concat_input_space_A * ...
                        opt_U_realizations <= concat_input_space_b + my_eps);
                    approx_stoch_u =sum(bin_u_orig)/n_particles;            
                    fprintf(['Input constraint satisfaction probability: ',...
                        '%1.3f | Acceptable: >%1.3f \n'], approx_stoch_u,...
                        1 - options.max_input_viol_prob);
                    fprintf(...
                        ['Undersampled probability (with %d particles): ',...
                        '%1.4f\nUnderapproximation to the original MILP (',...
                        'with %d particles): %1.4f\nAfter correction for',...
                        ' saturation: %1.4f\n'],...
                        n_kmeans, approx_voronoi_stoch_reach,...
                        n_particles, approx_stoch_reach_noc,approx_stoch_reach);
                end
            otherwise
        end
        kmeans_info.n_particles = n_particles;
        kmeans_info.n_kmeans = n_kmeans;
        kmeans_info.W_centroids = W_centroids;
        kmeans_info.W_realizations = W_realizations;
    end
end

function otherInputHandling(options)
    if ~(strcmpi(options.prob_str, 'term') &&...
            strcmpi(options.method_str, 'voronoi-affine'))
        throwAsCaller(SrtInvalidArgsError('Invalid options provided'));
    end
end
