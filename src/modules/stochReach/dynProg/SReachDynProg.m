function [prob_x, varargout] = SReachDynProg(prob_str, sys, x_inc, u_inc, ...
    safety_tube, varargin)
% Dynamic programming solution to stochastic reachability problems
% ============================================================================
%
% The function computes the probability of staying in a target tube defined
% on a particular state stace grid. SReachTools current REQUIRES the system to
% be LINEAR TIME-INVARIANT. The dynamic programming recursion can be found in 
%   
% S. Summers and J. Lygeros, "Verification of discrete time stochastic hybrid 
% systems: A stochastic reach-avoid decision problem," Automatica, vol. 46,
% no. 12, pp. 1951--1961, 2010.
%
% A trivial extension of this work to the case of time-varying safe set is
% implemented here.
%
% See also examples/doubleIntegratorDynamicProgramming.m.
%
% ============================================================================
%
% prob_x = SReachDynProg('term',sys,x_inc,u_inc,safety_tube)
% 
% Inputs:
% -------
%   prob_str    - String specifying the problem of interest. For each case, we
%                 compute the optimal value function that maps initial states
%                 to different maximal reach probabilities
%                     1. 'term' : Stay within the safety_tube
%   sys         - System description as a LtiSystem object
%   x_inc       - Scalar increment for all dimensions of the state space
%   u_inc       - Scalar increment for all dimensions of the input space
%   safety_tube - Safety tube of length N+1 where N is the time_horizon. It
%                 should have polyhedrons T_0, T_1, ...,T_N.
%
% Outputs:
% --------
%   prob_x      - Probability values at each grid point (M number of them) in
%                 grid_X (Mx1 array)
%   cell_xvec   - [Optional] Gridding along the particular dimension 
%                 (sys.state_dim x 1 cell array, with grid info along each
%                 dimension)
%   grid_x      - [Optional] Collection of grid points (Mx1 array)
%   mat_prob_x  - [Optional] M*(N+1) matrix of probability values corresponding
%                 to the "unrolled" value functions [V_0, V_1, ... V_N] where N
%                 is the time horizon. Note that prob_x = mat_prob_x(1,:)
%
% See also getDynProgLevelSets2D
%
% Notes:
% ------
% * REQUIRES:
%   - Gaussian-perturbed LtiSystem
%   - Input space is an axis-aligned HYPERCUBOID.
%   - State space is the smallest axis-aligned HYPERCUBOID that contains all the
%     sets in the target-tube
% * We impose uniform gridding across every dimension for the state and the
%   input.
% * WARNING: Dynamic programming suffers from the curse of dimensionality!
%   Using fine grids will increase the computation time.
% * SReachDynProg has a hidden `memoryusage` and `verbose` options. In future
%   versions, these will be handled via a `SReachDynProgOptions` struct.
%   - memoryusage governs the interplay between runtime and memory requirements
%     of dynamic programming
%        - memoryusage = 'high'
%            - Original behavior of SReachDynProg
%            - Compute the entire transition probability for all
%              (current_state, current_input, future_state) and then go
%              through the recursions. While this will lead to insanely fast
%              recursions, it will be memory intensive.
%        - memoryusage = 'low'
%            - Compute the entire transition probability for a given
%              current_state at every time step again and again. This will
%              lead to slower recursions, but it requires significantly
%              lesser memory.
%   - verbosity = {0,1} where verbose=0 implies quiet implementation and =1
%     provides feedback on progress of the dynamic programming
% ============================================================================
% 
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://sreachtools.github.io/license/

    % Input parsing
    valid_prob_str = {'term'}; %TODO: 'first'
    inpar = inputParser();
    inpar.addRequired('prob_str', @(x) any(validatestring(x,valid_prob_str)));
    inpar.addRequired('sys', @(x) validateattributes(x, ...
        {'LtiSystem'}, {'nonempty'})); %TODO: 'LtvSystem'
    inpar.addRequired('x_inc', @(x) validateattributes(x, {'numeric'}, ...
        {'scalar', '>', 0}));
    inpar.addRequired('u_inc', @(x) validateattributes(x, {'numeric'}, ...
        {'scalar', '>', 0}));
    inpar.addRequired('safety_tube',@(x) validateattributes(x,{'Tube'}, ...
        {'nonempty'}));

    try
        inpar.parse(prob_str, sys, x_inc, u_inc, safety_tube);
    catch err
        exc = SrtInvalidArgsError.withFunctionName();
        exc = exc.addCause(err);
        throwAsCaller(exc);
    end

    % 1. Ensure that the system is a Gaussian-perturbed LtiSystem
    % 2. Ensure that the state dim, input_dim <=4
    % 3. Input space is an axis-aligned hypercuboid
    % 4. safety_tube is appropriate
    % 5. optional arguments (in case of prob_str = 'first') is appropriate
    otherInputHandling(sys, safety_tube, prob_str, varargin);
    
    %% Hidden options for dynamic programming (TODO: create SReachDynOptions)
    % Switch between two memory usage options
    % low  - Recompute transition probability from every state to all
    %        possible input and next state combinations | Low memory usage
    % high - Compute all the transition probabilities at once => recursions
    %        are insanely fast | High memory usage
    memoryusage = 'low';   
    % Verbosity
    verbose = 0;
    % Update status every v_freq indices
    v_freq = 10;
    
    %% Compute the grid covering the target tube
    % n_targets is time_horizon + 1
    n_targets = length(safety_tube);    

    % Compute corners for gridding ==> Get corners of the largest target set in
    % safety_tube
    outerApproxVertices_target_sets = [];
    for itt = 1:n_targets
        outerApproxVertices_target_sets = [outerApproxVertices_target_sets;
            safety_tube(itt).outerApprox.V];
    end
    xmax = max(outerApproxVertices_target_sets);
    xmin = min(outerApproxVertices_target_sets);
    
    % Grid computation via allcomb
    if sys.state_dim == 1
        x1vec = xmin(1):x_inc:xmax(1);
        grid_x = x1vec';
        cell_xvec = {x1vec};
    elseif sys.state_dim == 2
        x1vec = xmin(1):x_inc:xmax(1);
        x2vec = xmin(2):x_inc:xmax(2);
        grid_x = allcomb(x1vec,x2vec);
        cell_xvec = {x1vec,x2vec};
    elseif sys.state_dim == 3
        x1vec = xmin(1):x_inc:xmax(1);
        x2vec = xmin(2):x_inc:xmax(2);
        x3vec = xmin(3):x_inc:xmax(3);
        grid_x = allcomb(x1vec,x2vec,x3vec);        
        cell_xvec = {x1vec,x2vec,x3vec};
    elseif sys.state_dim == 4
        x1vec = xmin(1):x_inc:xmax(1);
        x2vec = xmin(2):x_inc:xmax(2);
        x3vec = xmin(3):x_inc:xmax(3);
        x4vec = xmin(4):x_inc:xmax(4);
        grid_x = allcomb(x1vec,x2vec,x3vec,x4vec);        
        cell_xvec = {x1vec,x2vec,x3vec,x4vec};
    end
    % No. of grid points
    n_grid_x = length(grid_x);
    
    %% Input gridding
    umax = max(sys.input_space.V);
    umin = min(sys.input_space.V);
    if sys.input_dim == 1
        grid_u = allcomb(umin(1):u_inc:umax(1));
    elseif sys.input_dim == 2
        grid_u = allcomb(umin(1):u_inc:umax(1), ...
                         umin(2):u_inc:umax(2));
    elseif sys.input_dim == 3
        grid_u = allcomb(umin(1):u_inc:umax(1), ...
                         umin(2):u_inc:umax(2), ...
                         umin(3):u_inc:umax(3));
    elseif sys.input_dim == 4
        grid_u = allcomb(umin(1):u_inc:umax(1), ...
                         umin(2):u_inc:umax(2), ...
                         umin(3):u_inc:umax(3), ...
                         umin(4):u_inc:umax(4));
    end
    
    %% For trapezoid rule, we penalize 1/2 per dimension at the endpoints
    % Endpoint iff one of the dimensions is xmin or xmax
    % Use max to implement OR and sum to count how many active dimensions
    n_active_dims = max(sum(grid_x == xmin,2),sum(grid_x == xmax,2));
    % Compute scaling for the trapezoidal rule
    fraction_at_grid = 2.^(-n_active_dims);        
    % Create dx with appropriate scaling based on where it is located
    delta_x_grid = (x_inc^sys.state_dim).* fraction_at_grid;
    
    
    %% Initialize a matrix to store the value functions
    mat_prob_x = zeros(n_targets, n_grid_x); 
    
    %% How do we compute the transition probability
    switch memoryusage
        case 'high'
            % Compute transition probabilities for every x, u combination
            transition_prob_with_delta_all = computeTransProbWithDelta(sys, ...
                grid_x, grid_u, delta_x_grid, verbose, v_freq);
            % Return a n_grid_x * n_grid_u matrix corresponding current state 
            % grid_x(ix)
            transition_prob_with_delta = @(ix) ...
                transition_prob_with_delta_all{ix}';
        case 'low'
            % Compute the transition probability at every step
            transition_prob_with_delta = @(ix) ...
                computeTransProbWithDeltaAtX(sys, ix, grid_x, grid_u, ...
                    delta_x_grid)';
    end
    
    %% Implement dynamic programming 
    switch prob_str
        case 'term'
            terminal_indicator_x = safety_tube(n_targets).contains(grid_x');
            % fprintf('Set optimal value function at t=%d\n',n_targets-1);    
            mat_prob_x(n_targets,:) = terminal_indicator_x;
  
            for itt = n_targets - 1:-1:1
                % fprintf('Compute optimal value function at t=%d\n', itt - 1);
                % Obtain V_{t+1}
                old_prob_x = mat_prob_x(itt+1,:);
                % Check which of the grid points need to be iterated over
                current_indicator_x = safety_tube(itt).contains(grid_x');
                % Verbosity: Initial display
                if verbose 
                    switch memoryusage
                        case 'low'
                            n_ix = nnz(current_indicator_x);
                            fprintf(['Time t = %d | Computing V_t(x)....', ...
                                '%3.2f%%'], itt - 1, round(0/n_ix*100,2));
                        case 'high'
                            fprintf('Recursion at t = %d\n', itt - 1);
                    end
                end
                % Iterate over all these points and compute 
                % max_u \int_X V_{t+1}(x_{t+1}(u)) 
                %                       * transitionProb(x_{t+1},u)dx_{t+1}                    
                for ix = find(current_indicator_x==1)                    
                    mat_prob_x(itt,ix) = max(...
                        old_prob_x*transition_prob_with_delta(ix));
                    % Verbosity: continue_display
                    if verbose && strcmpi(memoryusage,'low') && ~mod(ix, v_freq)
                        percent_completed = round(ix/n_ix*100,2);
                        if percent_completed < 10
                            fprintf('\b\b\b\b\b%3.2f%%', percent_completed);
                        else
                            fprintf('\b\b\b\b\b\b%3.2f%%', percent_completed);
                        end                        
                    end                    
                end        
                if verbose && strcmpi(memoryusage,'low')
                    fprintf('\b\b\b\b\b\b%3.2f%%\n', 100);
                end
            end
        case 'first'    
            target_set = varargin{1};
            target_set_indicator_x = target_set.contains(grid_x');
            indx_grid_points_outside_target_set = ~double(target_set_indicator_x);
            
            % Initialize
            % fprintf('Set optimal value function at t=%d\n',n_targets-1);    
            mat_prob_x(n_targets,:) = target_set_indicator_x;

            for itt = n_targets - 1:-1:1
                % fprintf('Compute optimal value function at t=%d\n', itt - 1);
                % Obtain V_{t+1}
                old_prob_x = mat_prob_x(itt+1,:);
                % Set points that reached the target set as one
                mat_prob_x(itt,:) = target_set_indicator_x;
                % Iterate over the remaining safe points and compute
                % max_u \int_X V_{t+1}(x_{t+1}(u))transitionProb(x_{t+1},u)dx_{t+1}
                current_indicator_x = indx_grid_points_outside_target_set.* ...
                    safety_tube(itt).contains(grid_x');
                for ix = find(current_indicator_x==1)
                    mat_prob_x(itt,ix) = max(old_prob_x * ...
                        transition_prob_with_delta{ix}');
                end          
            end
    end
    prob_x = mat_prob_x(1,:);
    varargout{1} = cell_xvec;
    varargout{2} = grid_x;
    varargout{3} = mat_prob_x;
end

function transition_prob_with_delta = computeTransProbWithDelta(sys, grid_x, ...
    grid_u, delta_x_grid, verbose, v_freq)
    % Internal function to compute the transition probability scaled by the
    % Delta_x term for integration
    %
    % Returns the transition probability as a cell array of n_grid_u * n_grid_x 
    % matrices. Length of the cell array is n_grid_x (all possible current
    % states).
   

    n_grid_x = length(grid_x);
    
    % Define transition_prob as a cell array
    transition_prob_with_delta = cell(n_grid_x,1);               

    if verbose
        fprintf('Computing transition probability....%3.2f%%', ...
            round(0/n_grid_x*100,2));
    end
    for ix = 1:n_grid_x
        % Compute the transition probability at x(ix)
        transition_prob_with_delta{ix} = computeTransProbWithDeltaAtX(sys, ...
            ix, grid_x, grid_u, delta_x_grid);
        if verbose && ~mod(ix, v_freq)
            percent_completed = round(ix/n_grid_x*100,2);
            if percent_completed < 10
                fprintf('\b\b\b\b\b%3.2f%%', percent_completed);
            else
                fprintf('\b\b\b\b\b\b%3.2f%%', percent_completed);
            end                        
        end                    
    end
    if verbose == 1
        fprintf('\b\b\b\b\b\b%3.2f%%\n', 100);
    end
end


function transition_prob_with_delta_at_x = computeTransProbWithDeltaAtX(sys, ...
    ix, grid_x, grid_u, delta_x_grid)

    % Get transition probability for all points
    transition_prob_with_delta_at_x = zeros(length(grid_u), length(grid_x));
    
    % Compute the random vector (x_{k+1} - B u_k), given by A x_k + F w_k
    next_x_zi = sys.state_mat * grid_x(ix,:)' + sys.dist_mat * sys.dist;   
    
    for iu = 1:length(grid_u)
        % Compute the random vector x_{k+1} under a specific u_k
        next_x = next_x_zi + sys.input_mat * grid_u(iu,:);
        % Transition probability is the pdf times delta for integration
        transition_prob_with_delta_at_x(iu,:) = mvnpdf(grid_x, ...
            next_x.mean()', next_x.cov())' .* delta_x_grid';
    end
end

function otherInputHandling(sys, safety_tube, prob_str, optional_args)
    % 1. Ensure that the system is a Gaussian-perturbed LtiSystem
    % 2. Ensure that the state dim, input_dim <=4
    % 3. Input space is an axis-aligned hypercuboid
    % 4. safety_tube is appropriate
    % 5. optional arguments (in case of prob_str = 'first') is appropriate
    
    % Ensure the system is a Gaussian-perturbed system
    if ~(isa(sys.dist,'RandomVector') && strcmpi(sys.dist.type, 'Gaussian'))
        throwAsCaller(SrtInvalidArgsError(['Expected a Gaussian-perturbed', ...
            ' LTI System']));
    end
    
    % We can handle only 4-dim input space
    if sys.input_dim >4
        throwAsCaller(SrtInvalidArgsError('System can have at most 4 inputs'));
    end
    
    % We can handle only 4-dim state space
    if sys.state_dim >4
        throwAsCaller(SrtInvalidArgsError('System can have at most 4 states'));
    end
    
    % In order to go to hypercuboid, we will have to identify the edges
    % correctly. Currently a simple 2^-no_of_active_dims is used.
    % We will use MPT3's built in outer approximation code to confirm if this is
    % indeed the case
    if sys.input_space ~= sys.input_space.outerApprox
        throwAsCaller(SrtInvalidArgsError('Expected axis-aligned input space'));
    end
        
    % Safety tube
    if ~(isa(safety_tube, 'Tube') && safety_tube.tube(1).Dim == sys.state_dim) 
        throwAsCaller(SrtInvalidArgsError(['Expected a safety_tube of', ...
            ' dimension sys.state_dim']));
    end
    
    % Ensure that a target set alone was provided
    if strcmpi(prob_str, 'first')
        if isempty(optional_args)
            % Need one extra input
            throwAsCaller(SrtInvalidArgsError(...
                'Expected a target set (Polyhedron)'));
        elseif length(optional_args) >= 2
            % Only one additional argument
            throwAsCaller(SrtInvalidArgsError('Too many input arguments'));
        else
            target_set = optional_args{1};
            % Ensure target_set is a non-empty Polyhedron
            if ~(isa(target_set, 'Polyhedron') && ~target_set.isEmptySet()... 
                && target_set.Dim == sys.state_dim)
                err = SrtInvalidArgsError(['Expected a non-empty ', ...
                    'polyhedron of dimension sys.state_dim as target set']);
                throw(err);
            end
        end
    elseif strcmpi(prob_str, 'term') && ~isempty(optional_args)
        % No additional argument
        throwAsCaller(SrtInvalidArgsError('Too many input arguments'));
    end
end
