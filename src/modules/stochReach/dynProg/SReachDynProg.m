function [prob_x, varargout] = SReachDynProg(prob_str, sys, x_inc, u_inc,...
    safety_tube, varargin)
% Dynamic programming solution to stochastic reachability problems
% ============================================================================
%
% The function computes the probability of staying in a target tube defined
% on a particular state stace grid. The dynamic programming recursion can be
% found in 
%   
% S. Summers and J. Lygeros, "Verification of discrete time stochastic hybrid 
% systems: A stochastic reach-avoid decision problem," Automatica, vol. 46,
% no. 12, pp. 1951--1961, 2010.
%
% A trivial extension of this work to the case of time-varying safe set is
% implemented here.
%
% Usage: See example/doubleIntegratorDynamicProgramming.mlx
%
% ============================================================================
%
% prob_x=SReachDynProg('term',sys,x_inc,u_inc,safety_tube)
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
%                 should have polyhedrons T_0, T_1,...,T_N.
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
%                 to the "unrolled" value functions [V_0, V_1,... V_N] where N
%                 is the time horizon. Note that prob_x = mat_prob_x(1,:)
%
% See also getDynProgLevelSets2D
%
% Notes:
% ------
% * REQUIRES:
%   - Input space is an axis-aligned HYPERCUBOID.
%   - State space is the smallest axis-aligned HYPERCUBOID that contains all the
%     sets in the target-tube
% * We impose uniform gridding across every dimension.
% * WARNING: Dynamic programming suffers from the curse of dimensionality!
%   Using fine grids will increase the computation time.
% 
% ============================================================================
% 
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE

    % Input parsing
    valid_prob_str = {'term'}; %TODO: 'first'
    inpar = inputParser();
    inpar.addRequired('prob_str', @(x) any(validatestring(x,valid_prob_str)));
    inpar.addRequired('sys', @(x) validateattributes(x,...
        {'LtiSystem', 'LtvSystem'}, {'nonempty'}));
    inpar.addRequired('x_inc', @(x) validateattributes(x, {'numeric'},...
        {'scalar', '>', 0}));
    inpar.addRequired('u_inc', @(x) validateattributes(x, {'numeric'},...
        {'scalar', '>', 0}));
    inpar.addRequired('safety_tube',@(x) validateattributes(x,{'TargetTube'},...
        {'nonempty'}));

    try
        inpar.parse(prob_str, sys, x_inc, u_inc, safety_tube);
    catch err
        exc = SrtInvalidArgsError.withFunctionName();
        exc = exc.addCause(err);
        throwAsCaller(exc);
    end

    % 1. Ensure that the system is a Gaussian-perturbed LtiSystem
    % 2. Ensure that the state dim, input_dim <=3
    % 3. Input space is an axis-aligned hypercuboid
    % 4. safety_tube is appropriate
    % 5. optional arguments (in case of prob_str = 'first') is appropriate
    otherInputHandling(sys, safety_tube, prob_str, varargin);
    
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
    end
    
    n_grid_x = length(grid_x);
    
    %% Invoking trapezoid rule, we need to penalize 1/2 per dimension for the
    %% endpoints
    % max to implement OR and sum to count how many active dimensions
    n_active_dims = max(sum(grid_x == xmin,2),sum(grid_x == xmax,2));
    fraction_at_grid = 2.^(-n_active_dims);        
    delta_x_grid = (x_inc^sys.state_dim).* fraction_at_grid;
    
    % Input gridding
    umax = max(sys.input_space.V);
    umin = min(sys.input_space.V);
    if sys.input_dim == 1
        grid_u = allcomb(umin(1):u_inc:umax(1));
    elseif sys.input_dim == 2
        grid_u = allcomb(umin(1):u_inc:umax(1),...
                         umin(2):u_inc:umax(2));
    elseif sys.input_dim == 3
        grid_u = allcomb(umin(1):u_inc:umax(1),...
                         umin(2):u_inc:umax(2),...
                         umin(3):u_inc:umax(3));
    end
    
    % Initialize
    mat_prob_x = zeros(n_targets, n_grid_x); 
    
    % Compute transition probabilities
    transition_prob_with_delta = computeTransProbWithDelta(sys, grid_x,...
        grid_u, delta_x_grid);
    
    switch prob_str
        case 'term'
            terminal_indicator_x = safety_tube(n_targets).contains(grid_x');
%             fprintf('Set optimal value function at t=%d\n',n_targets-1);    
            mat_prob_x(n_targets,:) = terminal_indicator_x;
  
            for itt = n_targets - 1:-1:1
%                 fprintf('Compute optimal value function at t=%d\n', itt - 1);
                % Obtain V_{t+1}
                old_prob_x = mat_prob_x(itt+1,:);
                % Check which of the grid points need to be iterated over
                current_indicator_x = safety_tube(itt).contains(grid_x');
                % Iterate over all these points and compute 
                % max_u \int_X V_{t+1}(x_{t+1}(u))transitionProb(x_{t+1},u)dx_{t+1}
                for ix = find(current_indicator_x==1)
                    mat_prob_x(itt,ix) = max(...
                        old_prob_x*transition_prob_with_delta{ix}');
                end        
            end
        case 'first'    
            target_set = varargin{1};
            target_set_indicator_x = target_set.contains(grid_x');
            indx_grid_points_outside_target_set = ~double(target_set_indicator_x);
            
            % Initialize
%             fprintf('Set optimal value function at t=%d\n',n_targets-1);    
            mat_prob_x(n_targets,:) = target_set_indicator_x;

            for itt = n_targets - 1:-1:1
%                 fprintf('Compute optimal value function at t=%d\n', itt - 1);
                % Obtain V_{t+1}
                old_prob_x = mat_prob_x(itt+1,:);
                % Set points that reached the target set as one
                mat_prob_x(itt,:) = target_set_indicator_x;
                % Iterate over the remaining safe points and compute
                % max_u \int_X V_{t+1}(x_{t+1}(u))transitionProb(x_{t+1},u)dx_{t+1}
                current_indicator_x = indx_grid_points_outside_target_set.*...
                    safety_tube(itt).contains(grid_x');
                for ix = find(current_indicator_x==1)
                    mat_prob_x(itt,ix) = max(old_prob_x*transition_prob_with_delta{ix}');
                end          
            end
    end
    prob_x = mat_prob_x(1,:);
    varargout{1} = cell_xvec;
    varargout{2} = grid_x;
    varargout{3} = mat_prob_x;
end

function transition_prob_with_delta = computeTransProbWithDelta(sys, grid_x,...
    grid_u, delta_x_grid)
    % Internal function to compute the transition probability scaled by the
    % Delta_x term for integration

    verbose = 0;
    
    n_grid_x = length(grid_x);
    % Define transition_prob as a cell array
    transition_prob_with_delta = cell(n_grid_x,1);
    
    % Covariance matrix for x_k+1 is dist_mat * Sigma_dist * dist_mat'
    dist_cov = sys.dist_mat * sys.dist.parameters.covariance * sys.dist_mat';            

    if verbose == 1
        % For printing stuff --- Create fixed markers in the index space
        fprintf('Compute transition probability...000%%');
        if n_grid_x < 100
            no_of_splits = 10;
        else
            no_of_splits = 100;
        end
        print_marker = linspace(1,n_grid_x,no_of_splits+1);
        print_marker(end) = print_marker(end)-1;
        print_marker_indx = 1;
        print_marker_val = (print_marker(2)-print_marker(1))/n_grid_x*100;
    end
    
    for ix = 1:n_grid_x
        transition_prob_with_delta{ix} = zeros(length(grid_u), n_grid_x);
        for iu = 1:length(grid_u)
            % Compute the mean from the point of interest
            dist_mean = sys.state_mat * grid_x(ix,:)' +...
                sys.input_mat * grid_u(iu,:) +...
                sys.dist_mat * sys.dist.parameters.mean;
            % Transition probability is the pdf times delta for integration
            transition_prob_with_delta{ix}(iu,:) = mvnpdf(grid_x, dist_mean', dist_cov)'...
                .* delta_x_grid';
        end
        if verbose == 1 && ix > print_marker(print_marker_indx)
            val = (print_marker_indx-1) * print_marker_val;
            fprintf('\b\b\b\b%3d%%', round(val))
            print_marker_indx = print_marker_indx + 1;
        end
    end
    if verbose == 1
        fprintf('\b\b\b\b%3d%%', 100)
        fprintf('\n');
    end
end

function otherInputHandling(sys, safety_tube, prob_str, optional_args)

    % Ensure the system is a Gaussian-perturbed system
    if ~(isa(sys.dist,'RandomVector') && strcmpi(sys.dist.type, 'Gaussian'))
        throwAsCaller(SrtInvalidArgsError(['Expected a Gaussian-perturbed',...
            ' LTI System']));
    end
    
    % We can handle only 3-dim input space
    if sys.input_dim >=4
        throwAsCaller(SrtInvalidArgsError('System can have at most 3 inputs'));
    end
    
    % We can handle only 3-dim state space
    if sys.state_dim >=4
        throwAsCaller(SrtInvalidArgsError('System can have at most 3 states'));
    end
    
    % In order to go to hypercuboid, we will have to identify the edges
    % correctly. Currently a simple 2^-no_of_active_dims is used.
    % We will use MPT3's built in outer approximation code to confirm if this is
    % indeed the case
    if sys.input_space ~= sys.input_space.outerApprox
        throwAsCaller(SrtInvalidArgsError('Expected axis-aligned input space'));
    end
        
    % Safety tube
    if ~(isa(safety_tube, 'TargetTube') &&...
            safety_tube.tube(1).Dim == sys.state_dim) 
        throwAsCaller(SrtInvalidArgsError(['Expected a safety_tube of',...
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
                err=SrtInvalidArgsError(['Expected a non-empty polyhedron ',...
                    'of dimension sys.state_dim as target set']);
                throw(err);
            end
        end
    elseif strcmpi(prob_str, 'term') && ~isempty(optional_args)
        % No additional argument
        throwAsCaller(SrtInvalidArgsError('Too many input arguments'));
    end
end

%% Things for the first hitting time problem
%
% prob_x=SReachDynProg('first',sys,x_inc,u_inc,safety_tube,target_set)
%
%                     1. 'first' : Stay within the safety_tube and reach the
%                                  target set early if possible
%                     2. 'term' : Stay within the safety_tube
%
%   target_set  - [Required for 'first'] Target set that needs to be reached at
%                 some point within the horizon
