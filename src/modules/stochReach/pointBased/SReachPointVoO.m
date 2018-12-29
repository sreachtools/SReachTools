function [approx_stoch_reach, opt_input_vec, kmeans_info] = SReachPointVoO(...
    sys, initial_state, safety_tube, options)
% Solve the problem of stochastic reachability of a target tube (a lower bound
% on the maximal reach probability and an open-loop controller synthesis) using
% undersampled particle filter control
% =============================================================================
%
% SReachPointVoO implements a mixed-integer linear program-based approximation
% to the stochastic reachability of a target tube problem. To improve the
% tractability, we undersample the particles by utilizing the kmeans
% clustering/Voronoi partitioning of the disturbance samples. 
%
% H. Sartipizadeh, A. Vinod,  B. Acikmese, and M. Oishi, "Voronoi
% Partition-based Scenario Reduction for Fast Sampling-based Stochastic
% Reachability Computation of LTI Systems", In Proc. Amer. Ctrl. Conf., 2019
% 
% In contrast to `particle-open' approach, implemented in SReachPointPaO and
% described in
%
% K. Lesser, M. Oishi, and R. Erwin, "Stochastic reachability for control of
% spacecraft relative motion," in IEEE Conference on Decision and Control (CDC),
% 2013,
%
% SReachPointVoO permits a user-defined upper-bound on the overapproximation
% error.
%
%    High-level desc.   : Sample scenarios based on the additive noise and solve
%                         a mixed-integer linear program to make the maximum
%                         number of scenarios satisfy the reachability
%                         objective.  In addition, we use Voronoi partition to
%                         drastically improve the tractability while preserving
%                         the underapproximation quality
%    Approximation      : Overapproximation bounded above by a user-specified
%                         tolerance
%    Controller type    : Open-loop controller that satisfies the hard input
%                         bounds
%    Optimality         : Optimal (w.r.t scenarios drawn) open-loop controller
%                         for the underapproximation problem 
%
% =============================================================================
%
% [approx_stoch_reach, opt_input_vec] = SReachPointVoO(sys, initial_state,...
%   safety_tube, options)
%
% Inputs:
% -------
%   sys          - System description (LtvSystem/LtiSystem object)
%   initial_state- Initial state for which the maximal reach probability must be
%                  evaluated (A numeric vector of dimension sys.state_dim)
%   safety_tube  - Collection of (potentially time-varying) safe sets that
%                  define the safe states (Tube object)
%   options      - Collection of user-specified options for 'voronoi-open'
%                  (Matlab struct created using SReachPointOptions)
%
% Outputs:
% --------
%   approx_stoch_reach 
%               - An approximation of the stochastic reachability of a target
%                 tube problem computed using undersampled particle control
%                 approach using kmeans. In contrast to `particle-open'
%                 approach, this approximation permits a user-defined
%                 upper-bound on the overapproximation error.                   
%   opt_input_vec
%               - Open-loop controller: column vector of dimension
%                 (sys.input_dim*N) x 1
%   kmeans_info - A MATLAB struct containing the information about partitioning
%                 of GW space. The struct contains the following info:
%                  n_particles - Number of particles based off Hoeffding's
%                                inequality
%                  n_kmeans    - Number of bins for kmeans clustering
%                  GW_centroids- Centroids obtained from kmeans clustering
%                  GW_realizations
%                              - Realizations for the random vector GW
%
% See also SReachPoint.
%
% Notes:
% * This function requires CVX with Gurobi as the backend solver for optimizing
%   the resulting mixed-integer linear program.
% * This function requires kmeans function which is part of MATLAB's
%   Statistical and Machine Learning toolbox.
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

    % Ensure that system is stochastic
    if ~isa(sys.dist,'RandomVector')
        throwAsCaller(SrtInvalidArgsError('Expected a stochastic system'));
    end
    
    % Target tubes has polyhedra T_0, T_1, ..., T_{time_horizon}
    time_horizon = length(safety_tube)-1;

    approx_stoch_reach = -1;
    opt_input_vec = nan(sys.input_dim * time_horizon,1);

    % Requires Gurobi since we are solving a MILP
    [default_solver, solvers_cvx] = cvx_solver;
    if ~(contains(default_solver,'Gurobi') ||...
        any(contains(solvers_cvx,'Gurobi')))
        warning('SReachTools:setup',['SReachPointVoO returns a trivial ', ...
            'result since Gurobi (a MILP solver) was not setup.']);
    else
        % Ensure options is good
        otherInputHandling(options);
        
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

        % Step 1: Compute the number of particles needed to meet the
        % specified tolerance
        n_particles = ceil(...
            -log(options.failure_risk)/(2*options.max_overapprox_err^2));        
        if n_particles > 1e5
            throwAsCaller(SrtInvalidArgsError(...
                 sprintf(['Requested problem parameters required > 1e5 ',...
                    'particles (%d).'], n_particles)));
        end
        if options.n_kmeans > 100
            warning('SReachTools:runtime',sprintf(['Particle control with ',...
                'more than 100 samples may cause computational problems.\n',...
                'Going to analyze %d samples.'], options.n_kmeans));
        elseif options.n_kmeans > n_particles
            warning('SReachTools:runtime', ['Fewer particles needed than ',...
                'the number of kmean centroids specified. Using all ',...
                'particles']);
            options.n_kmeans = n_particles;
        end
        
        if options.verbose >= 1
            fprintf(['Required number of particles: %1.4e | Samples ',...
                'used: %3d\n'], n_particles, options.n_kmeans);
            fprintf('Creating random variable realizations....');
        end        
        % Compute the stochasticity of the concatenated disturbance random vec
        GW = G * concat(sys.dist, time_horizon);        
        % Create realizations of GW arranged columnwise
        GW_realizations = GW.getRealizations(n_particles);

        if options.verbose >= 1
            fprintf('Done\n');
        end

        % Implementation of Problem 3 in  Sartipizadeh ACC 2019 (submitted)
        % Step 2: Compute the Voronoi centers with prescribed number of
        % bins --- Use k-means clustering
        % Transposed input since kmeans expects each data point row-wise
        if options.verbose >= 1
            fprintf('Using k-means for undersampling....');
        end        
        [idx, GW_centroids_output] = kmeans(GW_realizations',...
            options.n_kmeans, 'MaxIter',1000);        
        GW_centroids = GW_centroids_output';
        if options.verbose >= 1
            fprintf('Done\n');
        end
        % Step 3: Compute the buffers associated with each of the Voronoi
        % centers
        buffers = zeros(n_lin_state, options.n_kmeans);
        voronoi_count = zeros(options.n_kmeans,1);
        for idx_indx = 1:options.n_kmeans
            relv_GW_realizations = (idx==idx_indx);
            voronoi_count(idx_indx) = nnz(relv_GW_realizations);
            GW_realizations_indx = GW_realizations(:, relv_GW_realizations);
            GW_centroid_indx = GW_centroids(:,idx_indx);
            % Buffer: max_k concat_safety_tube * (GW^(k) - GW_centroid)
            disp_from_centroid = concat_safety_tube_A * ...
                (GW_realizations_indx - repmat(GW_centroid_indx, 1,...
                    voronoi_count(idx_indx)));
            buffers(:, idx_indx) = max(disp_from_centroid,[],2);
        end
        % Solve the undersampled MILP and obtain a lower bound to the
        % original MILP
        % Step 4a: Solve the undersampled MILP
        if options.verbose >= 1
            fprintf('Setting up CVX problem....');
        end
        cvx_begin
            if options.verbose >= 2
                cvx_quiet false
            else
                cvx_quiet true
            end
            cvx_solver Gurobi
            variable X_realizations(sys.state_dim * time_horizon, options.n_kmeans);
            variable U_vector(sys.input_dim * time_horizon,1);
            variable z(1,options.n_kmeans) binary;
            maximize ((z*voronoi_count)/n_particles)
            subject to
                X_realizations == repmat(Z * initial_state + H * U_vector, ...
                    1, options.n_kmeans) + GW_centroids;
                concat_input_space_A * U_vector <= concat_input_space_b;
                concat_safety_tube_A * X_realizations + buffers <= repmat( ...
                    concat_safety_tube_b, 1, options.n_kmeans) + ...
                    options.bigM * repmat(1-z,n_lin_state,1);
            if options.verbose >= 1
                fprintf('Done\nParsing and solving the MILP....');
            end
        cvx_end       
        if options.verbose >= 1
            fprintf('Done\n');
        end
        %% Overwrite the solutions
        switch cvx_status
            case {'Solved','Inaccurate/Solved'}
                approx_voronoi_stoch_reach = cvx_optval;
                opt_input_vec = U_vector; 
                % Step 4b: Improve upon the estimate by a refined counting
                % Solve the mixed-integer linear program
                opt_X_realizations = repmat( ...
                    Z * initial_state + H * U_vector, 1, n_particles) + ...
                    GW_realizations;
                % All by default does columnwise (A particle succeeds only if
                % it clears all hyperplanes --- has a column of ones)
                z_orig = all(concat_safety_tube_A * opt_X_realizations <=...
                    concat_safety_tube_b);
                approx_stoch_reach = sum(z_orig)/n_particles -...
                    options.max_overapprox_err;            
                if options.verbose >= 1
                    fprintf(...
                        ['Undersampled probability (with %d particles): ',...
                        '%1.3f\nUnderapproximation to the original MILP ',...
                        '(with %d particles): %1.3f\n'],...
                        options.n_kmeans, approx_voronoi_stoch_reach,...
                        n_particles, approx_stoch_reach);
                end
            otherwise
        end
        kmeans_info.n_particles = n_particles;
        kmeans_info.n_kmeans = options.n_kmeans;
        kmeans_info.GW_centroids = GW_centroids;
        kmeans_info.GW_realizations = GW_realizations;
    end
end

function otherInputHandling(options)
    if ~(strcmpi(options.prob_str, 'term') &&...
            strcmpi(options.method_str, 'voronoi-open'))
        throwAsCaller(SrtInvalidArgsError('Invalid options provided'));
    end
end
