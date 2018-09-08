function [prob_x, cell_of_xvec, varargout] = getDynProgSolForTargetTube(sys, x_inc, u_inc,...
    target_tube)
% SReachTools/stochasticReachAvoid/getDynProgSolForTargetTube Get dynamic 
% programming grid probability for reachability of target tube
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
% Usage: See example doubleIntegratorDynmaicProgramming.m
%
% ============================================================================
%
% grid_prob = getDynProgSolForTargetTube(sys, ...
%     x_inc,...
%     u_inc,...
%     target_tube)
% 
% Inputs:
% -------
%   sys         - LtiSystem object
%   x_inc       - Scalar increment for all dimensions of the state space
%   u_inc       - Scalar increment for all dimensions of the input space
%   target_tube - Target tube of length N+1 where N is the time_horizon. It
%                 should have polyhedrons T_0, T_1,...,T_N.
%
% Outputs:
% --------
%   prob_x      - Mx1 Array of probability values at each grid point in
%                 grid_X
%   cell_of_xvec- sys.state_dim x 1 cell array each containing grid points
%                 along the particular dimension | Useful for
%                 getDynProgLevelSets
%   grid_x      - Mx1 Array of grid points obtained using
%                 allcomb(x1[,x2,x3]) where [ ] indicates the optional
%                 argument
%   mat_prob_x  - M*(N+1) matrix of probability values corresponding to the
%                 "unrolled" value functions [V_0, V_1,... V_N]. Note that
%                 prob_x = mat_prob_x(1,:)
%
% Notes:
% ------
% * REQUIRES:
%   - Input space is an axis-aligned HYPERCUBOID.
%   - State space is the smallest axis-aligned HYPERCUBOID fitting the
%     target-tube
% * For simplicity, we allow only uniform gridding across every dimension.
% * WARNING: Dynamic programming suffers from the curse of dimensionality!
%   Using fine grids will increase the computation time.
% 
% ============================================================================
% 
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE

    % check inputs
    validateattributes(sys, {'LtiSystem'}, {'nonempty'})
    validateattributes(target_tube, {'TargetTube'}, {'nonempty'});
    validateattributes(x_inc, {'numeric'}, {'scalar'});
    validateattributes(u_inc, {'numeric'}, {'scalar'});
    validateattributes(sys.dist, {'RandomVector', 'StochasticDisturbance'}, ...
        {'nonempty'});
    if ~strcmpi(sys.dist.type, 'Gaussian')
        throwAsCaller(SrtInvalidArgsError('Handles only Gaussian disturbances'));
    end
    
    % In order to go to hypercuboid, we will have to identify the edges
    % correctly. Currently a simple 2^-no_of_active_dims is used.
    if sys.input_space ~= sys.input_space.outerApprox
        throwAsCaller(SrtInvalidArgsError('Handles only axis-aligned input space'));
    end
    
    % n_targets is time_horizon + 1
    n_targets = length(target_tube);    

    % Compute corners for gridding ==> Get corners of the largest target set in
    % target_tube
    outerApproxVertices_target_sets = [];
    for itt = 1:n_targets
        outerApproxVertices_target_sets = [outerApproxVertices_target_sets;
            target_tube(itt).outerApprox.V];
    end
    xmax = max(outerApproxVertices_target_sets);
    xmin = min(outerApproxVertices_target_sets);

    if sys.state_dim == 1
        x1vec = xmin(1):x_inc:xmax(1);
        grid_x = x1vec';
        cell_of_xvec = {x1vec};
    elseif sys.state_dim == 2
        x1vec = xmin(1):x_inc:xmax(1);
        x2vec = xmin(2):x_inc:xmax(2);
        grid_x = allcomb(x1vec,x2vec);
        cell_of_xvec = {x1vec,x2vec};
    elseif sys.state_dim == 3
        x1vec = xmin(1):x_inc:xmax(1);
        x2vec = xmin(2):x_inc:xmax(2);
        x3vec = xmin(3):x_inc:xmax(3);
        grid_x = allcomb(x1vec,x2vec,x3vec);        
        cell_of_xvec = {x1vec,x2vec,x3vec};
    else
        throwAsCaller(SrtInvalidArgsError('System can have at most 3 states'));
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
    else
        throwAsCaller(SrtInvalidArgsError('System can have at most 3 inputs'));
    end
    
    fprintf('Set optimal value function at t=%d\n',n_targets-1);    
    terminal_indicator_x = target_tube(n_targets).contains(grid_x');

    % Initialize
    mat_prob_x = zeros(n_targets, n_grid_x); 
    mat_prob_x(n_targets,:) = terminal_indicator_x;
    
    transition_prob_with_delta = computeTransProbWithDelta(sys, grid_x, grid_u, delta_x_grid);
    for itt = n_targets - 1:-1:1
        fprintf('Compute optimal value function at t=%d\n', itt - 1);
        % Obtain V_{t+1}
        old_prob_x = mat_prob_x(itt+1,:);
        % Check which of the grid points need to be iterated over
        current_indicator_x = target_tube(itt).contains(grid_x');
        % Iterate over all these points and compute 
        % max_u \int_X V_{t+1}(x_{t+1}(u))transitionProb(x_{t+1},u)dx_{t+1}
        for ix = find(current_indicator_x==1)
            mat_prob_x(itt,ix) = max(old_prob_x*transition_prob_with_delta{ix}');
        end        
    end
    prob_x = mat_prob_x(1,:);
    varargout{1} = grid_x;
    varargout{2} = mat_prob_x;
end

function transition_prob_with_delta = computeTransProbWithDelta(sys, grid_x, grid_u, delta_x_grid)
    % Internal function to compute the transition probability scaled by the
    % Delta_x term for integration

    n_grid_x = length(grid_x);
    % Define transition_prob as a cell array
    transition_prob_with_delta = cell(n_grid_x,1);
    
    % Covariance matrix for x_k+1 is dist_mat * Sigma_dist * dist_mat'
    dist_cov = sys.dist_mat * sys.dist.parameters.covariance * sys.dist_mat';            
    
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
        if ix > print_marker(print_marker_indx)
            val = (print_marker_indx-1) * print_marker_val;
            fprintf('\b\b\b\b%3d%%', round(val))
            print_marker_indx = print_marker_indx + 1;
        end
    end
    fprintf('\b\b\b\b%3d%%', 100)
    fprintf('\n');
end
