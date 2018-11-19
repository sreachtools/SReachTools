function [opt_locations, separation] = spreadPointsOnUnitSphere(n_dim,...
    n_points, verbose)

    %% Initialize
    x_iter_unnorm = mvnrnd(zeros(n_points, n_dim), eye(n_dim))';
    norm_val = norms(x_iter_unnorm,2);
    x_iter = x_iter_unnorm./norm_val;
    
    %% Difference of convex approach
    continue_condition = 1;
    slack_tol = 1e-8;
    r_prev = Inf;
    sum_slack_prev = Inf;
    cost_tol = 1e-5;
    tau_iter = 1;
    scaling_tau = 1.1;
    tau_max = 1e4;
    max_iter = 20;
    iter_count = 1;
    
    %% Till the convergence condition is met
    while continue_condition
        slack_indx = 1;
        if verbose
            fprintf('%2d. Setting up the CVX problem...\n', iter_count);
        end
        cvx_begin quiet
            variable x(n_dim, n_points);
            variable r;
            variable slack_var(n_points * (n_points-1)/2, 1);
            minimize (-separation + tau_iter * sum(slack_var));
            subject to
                slack_var >= 0;
                separation >= 0;
                norm(x(:, 1)) <= 1;
                for pt_indx_1 = 2:n_points
                    norm(x(:, pt_indx_1)) <= 1;
                    if verbose
                        fprintf(' %2d |',pt_indx_1);
                        if abs(pt_indx_1-round(n_points/2))<1e-6
                            fprintf('\n');
                        end
                    end
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
                end
        if verbose
            fprintf('\nSolving the CVX problem...');
            cvx_end
            disp('done');
        else
            cvx_end
        end
        if any(contains({'Solved','Inaccurate/Solved'},cvx_status))
            cost_prev = -r_prev + tau_iter * sum_slack_prev;            
            if verbose
                fprintf(['Status: %s | Sum of slack: %1.3e (should be ~0) | ',...
                    'Change in opt cost: %1.3e\n'],...
                    cvx_status, sum(slack_var),...
                    abs(cvx_optval - cost_prev));
                violate_count = zeros(n_points,n_points);
                slack_indx = 1;
                for pt_indx_1 = 2:n_points
                    for pt_indx_2 = 1:pt_indx_1-1
                        %if norm(x(:, cup_indx_1) - x(:, cup_indx_2)) <=...
                            %2 * r_cupcake
                        if slack_var(slack_indx) >= 1e-6
                            violate_count(pt_indx_1, pt_indx_2) = 1;
                            violate_count(pt_indx_2, pt_indx_1) = 1;
                        end
                        slack_indx = slack_indx + 1;
                    end
                end                
            end            
            x_iter = x;
            % STOP if (slack small enough or slack converged) OR max iterations
            % CONTINUE if not of above with ORs replaced with AND
            continue_condition = ((sum(slack_var) > slack_tol) ||...
                (abs(cvx_optval - cost_prev) > cost_tol)) && ...
                (iter_count < max_iter);         
            r_prev = r;
            sum_slack_prev = sum(slack_var);
            iter_count = iter_count + 1;            
        else
            % Impossible to reach here since we are using slack variables
            keyboard;
        end
        tau_iter = min(tau_iter * scaling_tau, tau_max);    
        if verbose
            disp(' ');
        end
    end
    if iter_count > max_iter || (sum(slack_var) > slack_tol) ||...
            (abs(cvx_optval - cost_prev) > cost_tol)
        % Reached the maximum iteration but slack still not within
        % tolerance
        opt_locations = nan(2, n_points);
    else
        norm_val = norms(x,2);
        opt_locations = x./norm_val;
    end
end

% %% Check arrangment via difference of convex programming    
% clear;close all;clc;cvx_clear;
% 
% %% Check for a particular arrangement
% verbose = 1;
% n_dim = 3;
% n_points = 20;
% 
% [opt_locations, separation] = spread_points_on_sphere_dc(n_dim, n_points,...
%     verbose);
% 
% if n_dim == 2 || n_dim == 3
%     figure();
%     if n_dim == 2
%         quiver(zeros(1,n_points), zeros(1,n_points),...
%             opt_locations(1,:),opt_locations(2,:));
%     else
%         quiver3(zeros(1,n_points),zeros(1,n_points),zeros(1,n_points),...
%             opt_locations(1,:),opt_locations(2,:),opt_locations(3,:));
%     end
%     axis equal;
%     grid on;
%     box on;
% end
