function [opt_locations, separation] = spreadPointsOnUnitSphere(n_dim,...
    n_points, verbose)
% Computes a collection of n_points vectors of dimension n_dim that have a large
% minimum pairwise-separation
% ============================================================================
% 
% Given n_points, we spread n_points_first_quad = ceil(n_points - 2*n_dim) in
% the first quadrant using the following non-convex optimization problem,
% 
% maximize r
% subject to
%                  r >= 0,                                                   (1)
%  || x_i - x_j ||_2 >= r,   i,j\in{1,..,n_points_first_quad}, i<j           (2) 
%  || x_i - e_j ||_2 >= r,   i\in{1,..,n_points_first_quad}, j\in{1,...n_dim}(3)   
%        || x_i ||_2 <= 1,   i\in{1,..,n_points_first_quad}                  (4)  
%        || x_i ||_2 >= 0.8, i\in{1,..,n_points_first_quad}                  (5)   
%                x_i >= r/2                                                  (6)  
% where e_j refers to the standard vector (zeros with 1 at jth position). Here,
% (1) enforces positive separation, (2) and (3) enforces the smallest
% pairwise separation is above r (among each other and the standard vectors),
% (4) and (5) approximates || x_i ||_2 = 1, and (6) enforces the separation
% constraint is satisfied even among the reflections/rotations.
%
% This optimization problem is non-convex, and we solve it to a local optimality
% using difference-of-convex approach. Specifically, constraints (2), (3), and
% (5), which are reverse-convex, are tightened to their first-order Taylor
% series (under)approximation and the resulting linear constraints are enforced
% in their place. This method is discussed in:
%
%      J. D. Gleason, A. P. Vinod, and M. M. K. Oishi. 2018. Lagrangian 
%      Approximations for Stochastic Reachability of a Target Tube. 
%      online. (2018). https://arxiv.org/abs/1810.07118
%
% Next, we reflect/rotate these vectors in the first quadrant to occupy in all
% other quadrants. In the end, we tack on e_j and -e_j for each dimension j
% (hence, the - 2*n_dim).
%
% ============================================================================
% 
% [opt_locations, separation] = spreadPointsOnUnitSphere(n_dim,n_points,verbose)
%
% Inputs:
% -------
% n_dim     - Dimension of the unit sphere on which we wish to spread the points 
% n_points  - Number of points we wish to spread the points (Will be rounded up
%             to the smallest k such that 2^(n_dim) k + 2*n_dim >= n_points
% verbose   - Verbosity of this function
%
% Outputs:
% --------
% opt_locations - Unit vectors given as a n_dim x n_points
% separation    - The minimum pairwise-separation across the vectors
%
% Notes:
% ------
% * We enforce x_i^T x_i >= 0.8^2 instead of x_i^T x_i >= 1, so that the problem
%   converges faster.
% ============================================================================
%
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://sreachtools.github.io/license/
%
%

    %% Difference of convex approach
    continue_condition = 1;
    slack_tol = 1e-8;
    cost_tol = 1e-5;
    tau_iter = 1;
    scaling_tau = 1.1;
    tau_max = 1e4;
    max_iter = 20;
    iter_count = 1;
    
    n_points_first_quad = (n_points - 2*n_dim)/(2^n_dim);
    if n_points < 2*n_dim
        % Can't even allow standard vectors and their reflections? Unacceptable!
        throw(SrtInvalidArgsError('Expected n_points > 2*n_dim'));
    elseif mod(n_points - 2*n_dim,(2^n_dim)) > 0
        % Modify n_points such that it is 2^n_dim * k + 2*n_dim for some k> 0
        warning('SReachTools:runTime', sprintf(['Expected n_points = ',...
            '2^n_dim * k + 2*n_dim for some k> 0 | Got %d'], n_points));
        n_points_first_quad = ceil(n_points_first_quad);
        n_points = (2^n_dim) * n_points_first_quad + 2*n_dim;
    end
    if verbose
        fprintf('Spreading %d unit-length vectors in %d-dim space\n',...
            n_points, n_dim);
        fprintf('Analyzing %d unit-length vectors in first quadrant\n',...
            n_points_first_quad);
    end
    
    if n_points_first_quad > 0
        %% Initialize
        % Draw points from the multi-variate Gaussian and normalize it
        x_iter_unnorm = abs(mvnrnd(zeros(n_points_first_quad, n_dim),...
            10 * eye(n_dim))');
        norm_val = norms(x_iter_unnorm,2);
        x_iter = x_iter_unnorm./norm_val;
        % For difference of convex program, initialize the previous costs
        separation_prev = 0;
        sum_slack_prev = 0;

        %% Count the number of slack constraints (each category): 
        % Pairwise separation (2)
        n_pairwise_const = nchoosek(n_points_first_quad,2);
        % Pairwise separation from standard vectors (3)
        n_pairwise_const_plus_standard = n_pairwise_const + n_dim;
        % Norm reverse equality constraint (5)
        slack_count = n_pairwise_const_plus_standard + n_points_first_quad;

        %% Till the difference-of-convex convergence condition is met
        while continue_condition
            if verbose
                fprintf('%2d. Setting up the CVX problem...\n', iter_count);
            end
            cvx_begin quiet
                variable x(n_dim, n_points_first_quad);
                variable separation;
                variable slack_var(slack_count,1);
                minimize (-separation + tau_iter * sum(slack_var));
                subject to
                    slack_var >= 0;
                    separation >= 0;
                    slack_indx = 1;                 % Counter for slack variable
                    for pt_indx_1 = 1:n_points_first_quad
                        if verbose
                            fprintf(' %2d |',pt_indx_1);
                            if mod(pt_indx_1,10) < 1 &&... % == 0 check
                                    pt_indx_1 < n_points_first_quad
                                fprintf('\n');
                            end
                        end
                        % Constraint 2: Enforces the pairwise separation among 
                        % each other ||x_i - x_j||_2 >= r^2
                        for pt_indx_2 = 1:pt_indx_1-1
                            (x_iter(:, pt_indx_1) - x_iter(:, pt_indx_2))'*...
                               (x_iter(:, pt_indx_1) - x_iter(:, pt_indx_2)) + ...
                               2*[x_iter(:, pt_indx_1)-x_iter(:, pt_indx_2);
                                -(x_iter(:, pt_indx_1)-x_iter(:, pt_indx_2))]'*...
                               [x(:,pt_indx_1) - x_iter(:,pt_indx_1);
                                x(:,pt_indx_2) - x_iter(:,pt_indx_2)] +...
                               slack_var(slack_indx) >= separation^2;                        
                            slack_indx = slack_indx + 1;                    
                        end
                        % Constraint 3: Enforces the pairwise separation from 
                        % the standard vectors ||x_i - e_j||_2 >= r^2
                        for dim_indx = 1:n_dim
                            e_i_vector = zeros(n_dim,1);
                            e_i_vector(dim_indx) = 1;
                            (x_iter(:, pt_indx_1) - e_i_vector)'*...
                               (x_iter(:, pt_indx_1) - e_i_vector) + ...
                               2*(x_iter(:, pt_indx_1) - e_i_vector)'*...
                                 (x(:,pt_indx_1) - x_iter(:,pt_indx_1)) +...
                               slack_var(n_pairwise_const + dim_indx)...
                                >= separation^2;                        
                        end
                        % Constraint 4: Enforces the maximum norm constraint
                        % ||x_i||_2 <= 1
                        norm(x(:, pt_indx_1)) <= 1;
                        % Constraint 5: Enforces the minimum norm constraint
                        % ||x_i||_2 >= 0.8
                        x_iter(:, pt_indx_1)' * x_iter(:, pt_indx_1) +...
                            2 * x_iter(:,pt_indx_1)' *...
                            (x(:,pt_indx_1) - x_iter(:,pt_indx_1)) +...
                                slack_var(n_pairwise_const_plus_standard +...
                                    pt_indx_1) >= 0.8^2;                    
                        % Constraint 6: Enforces the separation constraint is
                        % satisfied even among the reflections/rotations
                        x >= separation/2;
                    end
            if verbose
                fprintf('\nSolving the CVX problem...');
                cvx_end
                disp('done');
            else
                cvx_end
            end
            switch cvx_status
                case {'Solved','Inaccurate/Solved'}
                    cost_prev = -separation_prev + tau_iter * sum_slack_prev;            
                    if verbose
                        fprintf(['Status: %s\nSum of slack: %1.3e ', ...
                            '(< %1.3e)\nChange in opt cost: %1.3e ', ...
                            '(< %1.3e)\n\n'],...
                            cvx_status, sum(slack_var), slack_tol, ...
                            abs(cvx_optval - cost_prev), cost_tol);
                    end            
                    x_iter = x;
                    % STOP if (slack small enough or slack converged) OR max 
                    % iterations
                    % CONTINUE if not of above with ORs replaced with AND
                    continue_condition = ((sum(slack_var) > slack_tol) ||...
                        (abs(cvx_optval - cost_prev) > cost_tol)) && ...
                        (iter_count < max_iter);         
                otherwise
                    % Impossible to reach here since we are using slack 
                    % variables
                    throw(SrtDevError(['Shouldn''t have reached here, ', ...
                        'since we are using slack variables']));
            end
            % Update the iteration values and other things for the next 
            % iteration
            iter_count = iter_count + 1;            
            tau_iter = min(tau_iter * scaling_tau, tau_max);    
            separation_prev = separation;
            sum_slack_prev = sum(slack_var);
        end
        if iter_count > max_iter || (sum(slack_var) > slack_tol) ||...
                (abs(cvx_optval - cost_prev) > cost_tol)
            % Reached the maximum iteration but slack still not within
            % tolerance
            opt_locations = nan(n_dim, n_points);
            warning('SReachTools:runTime',['Difference-of-convex program ',...
                'did not converge. Returning NaNs!']);
        else
            % Normalize the optimal locations
            norm_val = norms(x,2);
            opt_locations_first_quad = x./norm_val;
            % Reflect/rotate it to all the quadrants
            % sign_vectors has rows of [1 -1 combinations]
            sign_vectors = Polyhedron('lb',-ones(n_dim,1),'ub',ones(n_dim,1)).V;
            opt_locations = zeros(n_dim, n_points);
            opt_locations(:, 1:2*n_dim) = [eye(n_dim), -eye(n_dim)];
            for sign_indx = 1:size(sign_vectors,1)
                indx_locations = 2*n_dim + (sign_indx-1) * ...
                    n_points_first_quad+1:2*n_dim + ...
                        sign_indx*n_points_first_quad;
                opt_locations(:, indx_locations) = ...
                    diag(sign_vectors(sign_indx,:)) * opt_locations_first_quad;
            end     

        end
        if verbose
            disp('Completed spreading the vectors!');
        end
    else
        % Return the standard axis
        opt_locations(:, 1:2*n_dim) = [eye(n_dim), -eye(n_dim)];            
        separation = sqrt(2);
        if verbose
            disp('Skipped spacing vectors! Returned the standard axes!');
        end
    end
end
