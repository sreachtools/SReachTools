function varargout= generateMonteCarloSims(n_monte_carlo_sims, sys, ...
    initial_state, time_horizon, varargin)
% Generate Monte-Carlo simulations for (controlled/uncontrolled) LTI/LTV system 
% =============================================================================
% 
% generateMonteCarloSims produces a required number of trajectories,
% n_monte_carlo_sims, for a (affine-controlled/uncontrolled) LTI/LTV system sys
% with a deterministic/RandomVector initial_state for a given time_horizon. 
%
% If the system is controlled, then a causal disturbance-feedback affine
% controller may be specified (dist_feedback_gain, concat_input_vector). The
% controller will be saturated to the sys.input_space using projection. 
%
% For an open-loop controller, only an concat_input_vector may be specified and
% dist_feedback_gain is set to zero. 
%
% See also examples/forwardStochasticReachCWH.m, examples/cwhSReachPoint.m
%
% =============================================================================
% [concat_state_realization, concat_disturb_realizations, saturation_indx] = ...
%    generateMonteCarloSims(n_monte_carlo_sims, sys, initial_state, ...
%       time_horizon, optimal_input_vector, optimal_input_gain)
%
% Inputs:
% -------
%   n_monte_carlo_sims  - Number of Monte-Carlo simulation particles to be used 
%                         for estimation of the reach-avoid probability
%   sys                 - System description as a LtiSystem/LtvSystem object
%   initial_state       - Deterministic x_0
%   time_horizon        - Time horizon (N) of the stochastic reach-avoid problem
%   concat_input_vector - [Optional] Open-loop controller, a column vector of
%                         dimension (sys.input_dim*N) x 1 | Required only if 
%                         the system is controlled
%   dist_feedback_gain  - [Optional] Affine disturbance feedback gain for the
%                         concatenated disturbance vector, a matrix of dimension
%                         (sys.input_dim*N) x (sys.dist_dim*N) | Required only
%                         if the system is controlled, the controller is affine
%                         disturbance feedback, and the gain matrix must be
%                         lower block triangular (with zeros in its block
%                         diagonal elements) for causality | See Notes
%   srlcontrol          - [Optional but no concat_input_vector and empty
%                         dist_feedback_gain] SReachLagController object
%                         describing an admissible state feedback controller
%                         corresponding to the Lagrangian-based
%                         underapproximation to the stochastic reach set
%   verbose             - [Optional] Verbosity of this function when saturating
%                         affine disturbance feedback controllers
%
% Outputs:
% --------
%   concat_state_realization  - Matrix of concatenated state (column) vectors
%                               stacked columnwise. Each column has the state 
%                               trajectory [x_0; x_1; x_2; ...; x_N]
%   concat_disturb_realization- Matrix of concatenated disturbance (column) 
%                               vectors stacked columnwise. Each column has the 
%                               disturbance realization [w_0; w_1; ...; w_{N-1}]
%   saturation_indx           - [Available only for affine feedback] 
%                               Binary vector that indicates which realizations
%                               had their associated affine disturbance feedback
%                               controller saturated. Potentially non-zero only
%                               if the input_gain is non-empty
%   concat_input_realization  - [Available only for SReachLagController object] 
%                               Matrix of concatenated input (column) vectors
%                               stacked columnwise. Each column has the state 
%                               trajectory [x_0; x_1; x_2; ...; x_N]
%
% Notes:
% ------
% * Assumes IID disturbance for the LTI/LTV system. 
% * For controlled system, an open-loop controller NEEDS to be provided. The
%   optimal_input_vector should be a ((sys.input_dim) * time_horizon)-dim.
%   vector U = [u_0; u_1; ...; u_N] (column vector).
% * For uncontrolled system, the optimal_input_vector NEED NOT be provided
%   dist_feedback_gain must be lower
% * The disturbance feedback gain matrix must be lower block triangular, with
%   its block diagonal submatrices as zero matrices. This ensures that the
%   affine disturbance feedback controller's value at any point of time depends
%   only on the past disturbance values => causal controller. This function DOES
%   NOT check for this structure in the input gain.
% * Affine disturbance feedback controllers CAN NOT satisfy hard control bounds
%   when the disturbance is unbounded (like Gaussian). Therefore, we will
%   saturate the controller realization (associated with the disturbance
%   realization) via projection on to the concatenated input space.
%   Specifically, we solve the corresponding optimization problem for each
%   concatenated disturbance realization W
%   
%       minimize || U - (MW + D)||_2
%       subject to 
%           U \in \mathcal{U}^T
%
%   where U is the decision variable, \mathcal{U} is the input space, T is the
%   time horizon, M is the affine disturbance feedback gain, and D is the affine
%   disturbance feedback bias.
% * When using SReachLagController, verbosity may be specified in the following
%   way:
%       % Create a controller based on the underapproximation
%       srlcontrol = SReachLagController(sys, ... 
%           extra_info_under.bounded_dist_set, ...
%           extra_info_under.stoch_reach_tube);
%       % Generate Monte-Carlo simulations using the srlcontrol and
%       % generateMonteCarloSims
%       timer_mcarlo = tic;
%       [X,U,W] = generateMonteCarloSims(n_mcarlo_sims, sys, ...
%           initial_state, time_horizon, srlcontrol, [], ...
%           lagunder_options.verbose);
%
% ============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://sreachtools.github.io/license/
% 
%

    % Give progress reports in fractions defined by 1/disp_x when doing
    % saturation for affine controllers
    disp_x = 5;

    %% Input handling 
    % Ensure that n_monte_carlo_sims is a scalar and positive
    validateattributes(n_monte_carlo_sims, {'numeric'}, {'scalar',...
        'integer', '>', 0}, 'generateMonteCarloSims', 'n_monte_carlo_sims');
    % Ensure that sys is a LtiSystem/LtvSystem object
    validateattributes(sys, {'LtiSystem','LtvSystem'}, {'nonempty'},...
        'generateMonteCarloSims', 'sys');
    % Ensure initial state vector is valid or it is a random vector
    validateattributes(initial_state, {'RandomVector','numeric'},...
        {'nonempty'}, 'generateMonteCarloSims', 'initial_state');
    switch class(initial_state)
        case 'RandomVector'
            initial_state_realization = initial_state.getRealizations(1);
            if length(initial_state_realization) ~= sys.state_dim
                throwAsCaller(SrtInvalidArgsError(sprintf(['Expected a %d',...
                    '-dimensional random vector as initial_state'], ...
                    sys.state_dim)));
            end
        case 'numeric'
            if length(initial_state) ~= sys.state_dim
                throwAsCaller(SrtInvalidArgsError(sprintf(['Expected a %d',...
                    '-dimensional initial_state'], sys.state_dim)));
            end
    end
    % Ensure that time_horizon is a scalar and positive
    validateattributes(time_horizon, {'numeric'}, {'scalar',...
        'integer', '>', 0}, 'generateMonteCarloSims', 'time_horizon');
    
    % Input handling of optimal_input_vector and creation of
    % concat_optimal_input_vector
    concat_input_vector_repeated = zeros(time_horizon * sys.input_dim,...
        n_monte_carlo_sims);
    input_concat_dist_gain = zeros(time_horizon * sys.input_dim,...
        time_horizon * sys.dist_dim);
    if nargout > 3
        throwAsCaller(SrtInvalidArgsError('Too many output arguments'));
    elseif nargout == 3
        saturation_indx = zeros(n_monte_carlo_sims, 1);
    end

    verbose = 0;
    gain_present = 0;
    lag_control = 0;
    
    if sys.input_dim > 0 && ~isEmptySet(sys.input_space)
        if length(varargin) > 3
            throwAsCaller(SrtInvalidArgsError('Too many input arguments'));
        elseif ~isempty(varargin) && isa(varargin{1}, 'numeric')
            input_concat_vector = varargin{1};
            % Ensure optimal_input_vector is a valid open loop controller
            validateattributes(input_concat_vector, {'numeric'},...
                {'nonempty','vector','size',[time_horizon*sys.input_dim 1]},...
                'generateMonteCarloSims', 'input_concat_vector');
            % Repeat the same open-loop policy across all the simulations
            concat_input_vector_repeated = ...
                      repmat(input_concat_vector,1, n_monte_carlo_sims);
                  
            if length(varargin) >= 2
                if ~isempty(varargin{2})
                    gain_present = 1;
                    input_concat_dist_gain = varargin{2};
                    validateattributes(input_concat_dist_gain, {'numeric'},...
                        {'nonempty','size',[time_horizon * sys.input_dim,...
                                 time_horizon * sys.dist_dim]},...
                    'generateMonteCarloSims', 'input_concat_dist_gain');                
                else
                    input_concat_dist_gain= [];
                end
                % TODO: Check for lower block triangular structure
                if length(varargin) == 3
                    verbose = varargin{3};
                    validateattributes(verbose, {'numeric'},...
                        {'nonempty','scalar','integer','>=',0,'<=',1},...
                        'generateMonteCarloSims', 'verbose');
                end
            end
        elseif ~isempty(varargin) && isa(varargin{1}, 'SReachLagController')            
            lag_control = 1;
            % Input validation done at the elseif
            srlcontrol = varargin{1};
            if length(varargin) == 3 && isempty(varargin{2})
                verbose = varargin{3};
                validateattributes(verbose, {'numeric'},...
                    {'nonempty', 'scalar', 'integer', '>=', 0, '<=', 2}, ...
                    'generateMonteCarloSims', 'verbose');
            end
        else
            %warning('SReachTools:runTime','Setting input vectors to zero');
        end
    end

    W = sys.dist.concat(time_horizon);
    % Realization of the concatenated disturbance random vectors with each
    % realization stored columnwise
    if verbose >= 1
        fprintf('Getting %d realizations...', n_monte_carlo_sims);
    end
    concat_disturb_realizations = W.getRealizations(n_monte_carlo_sims);
    if verbose >= 1
        fprintf('Done\n');
    end    
        
    if lag_control == 1
        if verbose >= 1
            % disp('Rearranging W realizations');
            fprintf(['Pruning infeasible disturbance trajectories,\nwhen it', ... 
                ' violates at a particular time instant... %6d/%6d'], 0, ...
                n_monte_carlo_sims);
        end    
        
        % determine which realizations are "good", i.e. they satisfy the set
        % constraints
        good_rlz = all(reshape(srlcontrol.dist_set.contains( ...
            reshape(concat_disturb_realizations, sys.dist_dim, [])), [], ...
            n_monte_carlo_sims));
        
        % get the number of good trajectories and their indices in the
        % realization matrix
        n_traj = nnz(good_rlz);
        trajectory_indx = find(good_rlz);
        
        % initialization of state and input vectors
        concat_state_realizations = nan(sys.state_dim * (time_horizon + 1), ...
            n_traj);
        concat_input_realizations = zeros(sys.input_dim * time_horizon, n_traj);
        concat_state_realizations(1:sys.state_dim, :) = repmat(initial_state,...
            1, n_traj);
        if verbose >= 1
            fprintf(['\b\b\b\b\b\b\b\b\b\b\b\b\bDone.\nFound %d feasible ', ...
                'disturbance trajectories (~ %1.4f success probability)\n'], ...
                n_traj, n_traj/n_monte_carlo_sims);    
            fprintf('Rearranging realizations...\n');
            fprintf('Analyzing particles: %6d/%6d', 0, n_traj);
        end

        for t_indx = 1:n_traj
            W_realization = concat_disturb_realizations(:, ...
                trajectory_indx(t_indx));
            if verbose >= 1
                fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b%6d/%6d', t_indx, ...
                    n_traj);
            end

            for current_time = 1:time_horizon
                prev_time = current_time - 1;
                % Set the previous state
                prev_state = concat_state_realizations(...
                    prev_time*sys.state_dim+1:...
                    (prev_time+1)*sys.state_dim, t_indx);
                % Get the previous input based on the previous state and
                % the controller computable via the robust reach tube
                prev_action = srlcontrol.getInput(prev_state, prev_time);
                concat_input_realizations((prev_time)*sys.input_dim+1:...
                        (prev_time+1)*sys.input_dim) = prev_action;
                % Get the previous disturbance realization                
                prev_dist = W_realization(prev_time*sys.dist_dim+1:...
                        (prev_time+1)*sys.dist_dim);
                % Compute the current state
                concat_state_realizations((current_time)*sys.state_dim+1:...
                        (current_time+1)*sys.state_dim, t_indx) =...
                    sys.state_mat(prev_time) * prev_state +...
                        sys.input_mat(prev_time) * prev_action +...
                        sys.dist_mat(prev_time) * prev_dist;
            end
        end
        fprintf('\n');
        % Skip the initial state
        varargout{1} = concat_state_realizations(sys.state_dim+1:end,:);
        varargout{2} = concat_disturb_realizations;
        varargout{3} = concat_input_realizations;
    else
        % Compute the concatenated matrices for X
        [Z, H, G] = getConcatMats(sys, time_horizon);

        if isa(initial_state,'RandomVector')
            % Draw realizations from the initial state random vector
            concat_initial_state = initial_state.getRealizations( ...
                n_monte_carlo_sims);
        else
            % Concatenation of initial state vector
            concat_initial_state = repmat(initial_state,1, n_monte_carlo_sims);
        end

        % Realization of the random vector X (columnwise)
        % See @LtiSystem/getConcatMats for more info on the notation used
        if verbose
            fprintf(['Computing the reach probability associated with ', ...
                'the given controller via %1.2e Monte-Carlo simulation\n'], ...
                n_monte_carlo_sims);
        end

        if gain_present
            if verbose
                disp(['Affine disturbance feedback controller will be ', ...
                    'saturated to the input space via projection']);
            end
            concat_input_realizations = input_concat_dist_gain *...
                concat_disturb_realizations + concat_input_vector_repeated;
            % Check how many input realizations require saturation
            [concat_input_space_A, concat_input_space_b] = getConcatInputSpace(...
                sys, time_horizon);
            concat_input_space = Polyhedron('A', concat_input_space_A,...
                'b', concat_input_space_b);
            if verbose
                fprintf(['Using Polyhedron/contains to identify ', ...
                    'realizations that require saturation...']);
            end
            saturation_indx =...
                concat_input_space.contains(concat_input_realizations);        
            proj_req = find(saturation_indx == 0);
            if verbose
                disp('Done');
                fprintf('Input constraint violation probability: %1.4f\n',...
                    length(proj_req)/n_monte_carlo_sims);
                fprintf(['We need to saturate %d realizations. We will '...
                    'provide progress in %d quantiles.\n'], length(proj_req),...
                    disp_x);
            end

            % Saturate these realizations
            realization_counter = 0;            
            realization_frac_disp = round(length(proj_req)/disp_x); %disp every x
            for realization_indx = proj_req
                % Affine disturbance feedback controller associated with
                % concat_disturb_realization (U = MW + D)
                realization_counter = realization_counter + 1;            
                if verbose && ...
                   (mod(realization_counter, realization_frac_disp) == 0)

                    fprintf('Completed saturating %5d/%5d input realizations\n',...
                       realization_counter, length(proj_req));
                end
                % Saturate the resulting inputs
                cvx_begin quiet
                    variable concat_saturated_input(sys.input_dim * time_horizon, 1);

                    minimize (norm(concat_saturated_input -...
                        concat_input_realizations(:, realization_indx)))
                    subject to
                        concat_input_space_A * concat_saturated_input <=...
                            concat_input_space_b;                    
                cvx_end
                switch cvx_status
                    case {'Solved','Inaccurate/Solved'}
                        % Input realization saturated via projection
                        concat_input_realizations(:, realization_indx) =...
                            concat_saturated_input;
                        if cvx_optval > 1e-8
                            saturation_indx(:, realization_indx) = 1;
                        end
                    otherwise
                        throw(SrtDevError(['Saturation via projection ', ...
                            'failed. This might be due to numerical issues.']));
                end
            end
        else
            concat_input_realizations = concat_input_vector_repeated;
        end
        concat_state_realization = [concat_initial_state;
            Z * concat_initial_state + H * concat_input_realizations + ...
              G * concat_disturb_realizations];
        varargout{1} = concat_state_realization;
        varargout{2} = concat_disturb_realizations;
        varargout{3} = concat_input_realizations;
        if nargout > 3
            varargout{4} = saturation_indx;
        end
    end
end
