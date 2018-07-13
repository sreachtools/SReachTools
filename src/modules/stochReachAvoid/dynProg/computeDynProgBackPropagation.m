function grid_prob = computeDynProgBackPropagation(sys, ...
    state_grid, input_grid, grid_prob, initial_set)
% SReachTools/stochasticReachAvoid/computeDynProgBackPropagation Compute the
% dynamic programming back propagation
% ============================================================================
%
% The function computes the one-step back propagation for the dynamic 
% programming recursion. See
% 
% S. Summers and J. Lygeros, "Verification of discrete time stochastic hybrid 
% systems: A stochastic reach-avoid decision problem," Automatica, vol. 46,
% no. 12, pp. 1951--1961, 2010.
%
% Usage: See getDynProgSolForTargetTube
% 
% ============================================================================
%
% grid_prob = computeDynProgBackPropagation(sys, ...
%     state_grid, input_grid, grid_prob, initial_set)
% 
% Inputs:
% -------
%   sys         - LtiSystem object
%   state_grid  - SpaceGrid object
%   input_grid  - InputGrid object
%   grid_prob   - Nx1 Array of probability values, where N is equivalent
%                 to size(state_grid, 1)
%   initial_set - Polyhedron object
%
% Outputs:
% --------
%   grid_prob - Nx1 Array of probability values, where N is equivalent
%                      to size(state_grid, 1)
%
% See also getDynProgSolForTargetTube
%
% Notes:
% ------
% * Currently this back propagation, and subsequently the entire dynamic 
%   programming recursion, only works for Gaussian disturbances.
%
% ============================================================================
%
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 

    % save the original grid probability
    old_grid_prob = grid_prob;
    n_state_grid_points = size(state_grid.grid, 1);
    n_input_grid_points = size(input_grid.grid, 1);
    
    % for low dimensional systems we can speed up computation
    if state_grid.dim == 2
        ext_grid = state_grid.getExternalGrid();
        use_fast_gaussian = true;
    else
        use_fast_gaussian = false;
    end
    
    for ix = 1:n_state_grid_points
        state_vec = state_grid.grid(ix, :)';
        if ~initial_set.contains(state_vec)
            % since not in the initial set can immediately set probability to 
            % zero
            grid_prob(ix) = 0;
        else
            input_prob = zeros(n_input_grid_points, 1);
            for iu = 1:n_input_grid_points
                if use_fast_gaussian
                    % 2-d fase gaussian method
                    transition_prob = fastGaussianProbFor2d(sys, ext_grid, ...
                        state_vec, input_grid.grid(iu, :)');
                else
                    % normal gaussian method
                    transition_prob = computeGaussianProbabForInputAndState(...
                        sys, ...
                        state_grid, ...
                        state_vec, ...
                        input_grid.grid(iu, :)');
                end

                % probability for each input option
                input_prob(iu) = transition_prob' * old_grid_prob;
            end 

            % get the highest probability
            grid_prob(ix) = max(input_prob);
        end
    end

end

function probability = fastGaussianProbFor2d(sys, ext_grid, state_vec, ...
    input_vec)
% SReachTools/stochasticReachAvoid/fastGaussianProbFor2d
% ============================================================================
%
% The function computes the transition probabilites for a 2-d system with
% a Gaussian disturbance, leveraging external grids and diff operations
% 
% Usage: Nested function
%
% ============================================================================
%
% probability = fastGaussianProbFor2d(sys, ext_grid, state_vec, ...
%     input_vec)
% 
% Inputs:
% -------
%   sys       - LtiSystem object
%   ext_grid  - SpaceGrid object
%   stat_vec  - State vector of current point
%   input_vec - Current input vector
%
% Outputs:
% --------
%   probability - Probability of transisioning from state_vec -> point on the 
%                 original state grid
%
% Notes:
%   - Function is exclusively used within 
%     computeDynProgBackPropagation.
%
% ============================================================================
% 
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 

    % compute the mean of the Gaussian through the state transition
    mu = sys.state_mat * state_vec + ...
        sys.input_mat * input_vec + ...
        sys.dist_mat * sys.dist.parameters.mean;
    
    % variance after state transition
    sigma = sys.dist_mat * sys.dist.parameters.covariance * ...
        sys.dist_mat';
    
    % compute the probability at points on the original state grid
    probability = mvncdf(ext_grid.grid, mu', sigma);
    probability = reshape(probability, ext_grid.n_points);
    probability = diff(diff(probability,1,2));
    probability = reshape(probability,[],1);
end

function probability = computeGaussianProbabForInputAndState(sys, ...
    state_grid, state_vec, input_vec)
% SReachTools/computeDynProgBackPropagation/computeGaussianProbabForInputAndState
% Comput gaussian transition probability
% ============================================================================
%
% The function computes the transition probabilites for a 2-d system with
% a Gaussian disturbance, leveraging external grids and diff operations
%
% Usage: Nested function
% 
% ============================================================================
%
% probability = computeGaussianProbabForInputAndState(sys, ...
%     state_grid, state_vec, input_vec)
% 
% Inputs:
% -------
%   sys         - LtiSystem object
%   state_grid  - SpaceGrid object
%   stat_vec    - State vector of current point
%   input_vec   - Current input vector
%
% Outputs:
% --------
%   probability - Probability of transisioning from state_vec -> point on the 
%                 original state grid
%
% Notes:
%   - Function is exclusively used within 
%     computeDynProgBackPropagation.
% 
% ============================================================================
%
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 

    % mean after state transition
    mu = sys.state_mat * state_vec + ...
        sys.input_mat * input_vec + ...
        sys.dist_mat * sys.dist.parameters.mean;
    
    % covariance after state transition
    sigma = sys.dist_mat * sys.dist.parameters.covariance * ...
        sys.dist_mat';
    
    % initialize and comput probabilty on the state grid
    probability = zeros(size(state_grid.grid, 1), 1);
    for i = 1:length(probability)
        box = SimpleBox(point, dx);
        probability(i) = box.computeGaussianProbability(...
            mvncdf(box.vertices, mu', sigma));
    end
end

function p = computeProbabilityAtGridPoint(point, dx, mu, sigma)
% SReachTools/computeDynProgBackPropagation/computeProbabilityAtGridPoint
% Compute the probability at a grid point
% ============================================================================
%
% The function computes the transition probabilites for a 2-d system with
% a Gaussian disturbance, leveraging external grids and diff operations
% 
% Usage: Nested function
% 
% ============================================================================
%
% p = computeProbabilityAtGridPoint(point, dx, mu, sigma)
% 
% Inputs:
% -------
%   point - Current State Space Grid point
%   dx    - Grid spacing
%   mu    - Mean value of Gaussian after state transition
%   sigma - Covariance of Gaussian after state transition
%
% Outputs:
% --------
%   p - Transition probability to point
%
% Notes:
%   - Function is exclusively used within 
%     computeDynProgBackPropagation.
% 
% ============================================================================
% 
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 

    box = SimpleBox(point, dx);
    p = box.computeGaussianProbability(mvncdf(box.vertices, mu', sigma));
end

function new_p = gaussianProbabilityDifference(old_p, ...
    old_grid_points)
% SReachTools/computeDynProgBackPropagation/gaussianProbabilityDifference
% Compute the probability from gaussian using box and diffs
% ============================================================================
%
% The function takes a bounding box and computes, by using differences, the
% probability of the interior point
% 
% Usage: Nested function
% 
% ============================================================================
%
% new_p = gaussianProbabilityDifference(old_p, ...
%     old_grid_points)
% 
% Inputs:
% -------
%   old_p           - Probabilities at grid points
%   old_grid_points - Grid points
%
% Outputs:
%   new_p - Probability at gird points after differencing
%
% Notes:
%   - Function is exclusively used within 
%     computeDynProgBackPropagation.
%   - Need to re-examing the names of inputs and outputs to provide more clarity
% 
% ============================================================================
% 
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE


    if size(old_grid_points, 1) == 4
        new_p = diff(diff(reshape(old_p, [2, 2]), 1, 2));
    else
        n_dims = size(old_grid_points, 2);
        n_old_grid_points = size(old_grid_points, 1);
        new_grid_points = zeros(2^(n_dims-1), n_dims);
        n_new_grid_points = size(new_grid_points, 1);
        new_p = zeros(n_new_grid_points, 1);

        for i = 2:2:n_old_grid_points
            new_p = old_p(i-1) - old_p(i);
            new_grid_points = old_grid_points(i, 1:n_dims-1);
        end

        if size(new_p, 1) > 1
            new_p = gaussianProbabilityDifference(new_p, ...
                new_grid_points);
        end
    end
end

function box_points = getBoxPointsFromGridPoint(point, dx)
% SReachTools/computeDynProgBackPropagation/getBoxPointsFromGridPoint
% Get bounding box around a grid point
% ============================================================================
%
% This function computes points that create a bounding box around a given
% point with given spacing in each direction
%
% Usage: Nested function
% 
% ============================================================================
%
% box_points = getBoxPointsFromGridPoint(point, dx)
% 
% Inputs:
% -------
%   point - Current SpaceGrid grid point
%   dx    - State SpaceGrid deltas
%
% Outputs:
% --------
%   box_points - Array of grid points that creates a bounding box around 'point'
%                with spacings dx
%
% Notes:
%   - Function is exclusively used within 
%     computeDynProgBackPropagation.
% 
% ============================================================================
% 
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 

    if length(dx) == 2
        box_points = get2dBoxPointsFromGridPoint(point, dx);
    elseif length(dx) == 3
        box_points = get3dBoxPointsFromGridPoint(point, dx);
    else
        box_points = zeros(2^length(dx), size(point, 2));
        ones_vec   = [1, -1];
        ind_vec    = ones(1, length(dx));
        shift_vec  = zeros(1, length(dx));
        for i = 1:size(box_points, 1)
            for j = 1:length(dx)
                shift_vec(j) = ones_vec(ind_vec(j));
            end
            box_points(i, :) = point + shift_vec .* dx;

            ind_vec(end) = ind_vec(end) + 1;
            for j = length(dx):-1:1
                if ind_vec(j) > 2
                    if j == 1
                        break;
                    else
                        ind_vec(j) = 1;
                        ind_vec(j-1) = ind_vec(j-1) + 1;
                    end
                else
                    break;
                end
            end
        end
    end
    

end

function box_points = get2dBoxPointsFromGridPoint(point, dx)
% SReachTools/computeDynProgBackPropagation/get2dBoxPointsFromGridPoint
% Get 2d bounding box verties from grid point
% =============================================================================
%
% This function computes points that create a bounding box around a given
% point with given spacing in each direction for 2-dimensional systems. The
% small system can be hard-coded for increased speed.
%
% Usage: Nested function
% 
% =============================================================================
%
% box_points = get2dBoxPointsFromGridPoint(point, dx)
% 
% Inputs:
% -------
%   point - Current SpaceGrid grid point
%   dx    - State SpaceGrid deltas
%
% Outputs:
% --------
%   box_points - Array of grid points that creates a bounding box around 'point'
%                with spacings dx
%
% Notes:
%   - Function is exclusively used within 
%     computeDynProgBackPropagation.
% 
% =============================================================================
% 
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
%


    box_points = zeros(4, length(dx));
    box_points(1, :) = point + [1, 1] .* dx;
    box_points(2, :) = point + [1, -1] .* dx;
    box_points(3, :) = point + [-1, 1] .* dx;
    box_points(4, :) = point + [-1, -1] .* dx;
end

function box_points = get3dBoxPointsFromGridPoint(point, dx)
% SReachTools/computeDynProgBackPropagation/get3dBoxPointsFromGridPoint
% Get bounding box vertiices for 3-d box at grid point
% ============================================================================
%
% This function computes points that create a bounding box around a given
% point with given spacing in each direction for 3-dimensional systems. The
% small system can be hard-coded for increased speed.
%
% Usage: Nested function
% 
% ============================================================================
%
% box_points = get3dBoxPointsFromGridPoint(point, dx)
% 
% Inputs:
% -------
%   point - Current SpaceGrid grid point
%   dx    - State SpaceGrid deltas
%
% Outputs:
% --------
%   box_points - Array of grid points that creates a bounding box around 'point'
%                with spacings dx
%
% Notes:
%   - Function is exclusively used within 
%     computeDynProgBackPropagation.
% 
% ============================================================================
% 
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
%

    box_points = zeros(8, length(dx));
    box_points(1, :) = point + [1, 1, 1] .* dx;
    box_points(2, :) = point + [1, 1, -1] .* dx;
    box_points(3, :) = point + [1, -1, 1] .* dx;
    box_points(4, :) = point + [1, -1, -1] .* dx;
    box_points(5, :) = point + [-1, 1, 1] .* dx;
    box_points(6, :) = point + [-1, 1, -1] .* dx;
    box_points(7, :) = point + [-1, -1, 1] .* dx;
    box_points(8, :) = point + [-1, -1, -1] .* dx;
end
