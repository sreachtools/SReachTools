clear all
clc
close all

%% Monte-Carlo simulation parameters
n_mcarlo_sims = 1e5;
n_sims_to_plot = 5;

%% CWH system params
umax=0.1;
mean_disturbance = zeros(4,1);
covariance_disturbance = diag([1e-4, 1e-4, 5e-8, 5e-8]);
% Define the CWH (planar) dynamics of the deputy spacecraft relative to the chief spacecraft as a LtiSystem object
sys = getCwhLtiSystem(4,...
                      Polyhedron('lb', -umax*ones(2,1),...
                                 'ub',  umax*ones(2,1)),...
                      StochasticDisturbance('Gaussian',...
                                            mean_disturbance,...
                                            covariance_disturbance));
time_horizon=5;                                             % Stay within a line of sight cone for 4 time steps and 
                                                            % reach the target at t=5% Safe Set --- LoS cone
%% Safe set definition --- LoS cone |x|<=y and y\in[0,ymax] and |vx|<=vxmax and |vy|<=vymax
ymax=2;
vxmax=0.5;
vymax=0.5;
A_safe_set = [1, 1, 0, 0;           
             -1, 1, 0, 0; 
              0, -1, 0, 0;
              0, 0, 1,0;
              0, 0,-1,0;
              0, 0, 0,1;
              0, 0, 0,-1];
b_safe_set = [0;
              0;
              ymax;
              vxmax;
              vxmax;
              vymax;
              vymax];
safe_set = Polyhedron(A_safe_set, b_safe_set);
%% Target set --- Box [-0.1,0.1]x[-0.1,0]x[-0.01,0.01]x[-0.01,0.01]
target_set = Polyhedron('lb', [-0.1; -0.1; -0.01; -0.01],...
                        'ub', [0.1; 0; 0.01; 0.01]);
target_tube = TargetTube('reach-avoid',safe_set, target_set, time_horizon);                    
%% Problem 1 and 2: Verification and controller synthesis from a given initial state
% We will first specify the initial state and parameters for the MATLAB's Global 
% Optimization Toolbox |patternsearch|.
%%
initial_state = [-0.75;         % Initial x relative position
                 -0.75;         % Initial y relative position
                 0;             % Initial x relative velocity
                 0];            % Initial y relative velocity
slice_at_vx_vy = initial_state(3:4);             
%% Parameters for MATLAB's Global Optimization Toolbox patternsearch
desired_accuracy = 1e-3;        % Decrease for a more accurate lower 
                                % bound at the cost of higher 
                                % computation time
PSoptions = psoptimset('Display','off');

%% Generate matrices for optimal mean trajectory generation
[H_matrix, mean_X_sans_input, ~] =...
       getHmatMeanCovForXSansInput(sys, initial_state, time_horizon);
%% 
% Next, using SReachTools, we will compute an <https://doi.org/10.1109/LCSYS.2017.2716364 
% optimal open-controller and the associated reach-avoid probability>. This function 
% takes about few minutes to run.
timer_ft = tic;
disp('Genz-algorithm+patternsearch technique');
[lb_stoch_reach_avoid_ft, optimal_input_vector_ft] = ...
                         getLowerBoundStochReachAvoid(sys,...
                                                      initial_state,...
                                                      target_tube,...
                                                      safe_set,...
                                                      'genzps',...
                                                      [],...
                                                      desired_accuracy,...
                                                      PSoptions);  
elapsed_time_ft = toc(timer_ft);
% Monte-Carlo validation
% This function returns the concatenated state vector stacked columnwise
concat_state_realization_ft = generateMonteCarloSims(...
                               n_mcarlo_sims,...
                               sys,...
                               initial_state,...
                               time_horizon,...
                               optimal_input_vector_ft);
% Check if the location is within the target_set or not
mcarlo_result_ft = target_tube.contains(concat_state_realization_ft);
% Optimal mean trajectory generation                         
optimal_mean_X_ft = mean_X_sans_input + H_matrix * optimal_input_vector_ft;
optimal_mean_trajectory_ft=reshape(optimal_mean_X_ft,sys.state_dim,[]);                                              
    
%% CC (Iterative risk allocation)     
timer_cc_iter = tic;
disp('Iterative risk allocation technique');
[lb_stoch_reach_avoid_cc_iter, optimal_input_vector_cc_iter] =...
                         getLowerBoundStochReachAvoid(sys,...
                                                      initial_state,...
                                                      target_tube,...
                                                      safe_set,...
                                                      'ccciter',...
                                                      desired_accuracy);  
elapsed_time_cc_iter = toc(timer_cc_iter);
% This function returns the concatenated state vector stacked columnwise
concat_state_realization_cc_iter = generateMonteCarloSims(...
                               n_mcarlo_sims,...
                               sys,...
                               initial_state,...
                               time_horizon,...
                               optimal_input_vector_cc_iter);
% Check if the location is within the target_set or not
mcarlo_result_cc_iter = target_tube.contains(concat_state_realization_cc_iter);
% Optimal mean trajectory generation                         
optimal_mean_X_cc_iter = mean_X_sans_input + H_matrix * optimal_input_vector_cc_iter;
optimal_mean_trajectory_cc_iter=reshape(optimal_mean_X_cc_iter,sys.state_dim,[]);

%% CC (Linear program approach)     
timer_cc_pwl = tic;
disp('Piecewise-linear single shot risk allocation technique');
[lb_stoch_reach_avoid_cc_pwl, optimal_input_vector_cc_pwl] =...
                         getLowerBoundStochReachAvoid(sys,...
                                                      initial_state,...
                                                      target_tube,...
                                                      safe_set,...
                                                      'cccpwl',...
                                                      desired_accuracy);  
elapsed_time_cc_pwl = toc(timer_cc_pwl);
% This function returns the concatenated state vector stacked columnwise
concat_state_realization_cc_pwl = generateMonteCarloSims(...
                               n_mcarlo_sims,...
                               sys,...
                               initial_state,...
                               time_horizon,...
                               optimal_input_vector_cc_pwl);
% Check if the location is within the target_set or not
mcarlo_result_cc_pwl = target_tube.contains(concat_state_realization_cc_pwl);
% Optimal mean trajectory generation                         
optimal_mean_X_cc_pwl = mean_X_sans_input + H_matrix * optimal_input_vector_cc_pwl;
optimal_mean_trajectory_cc_pwl=reshape(optimal_mean_X_cc_pwl,sys.state_dim,[]);

%% Plotting
figure();
box on;
hold on;
plot(safe_set.slice([3,4], slice_at_vx_vy), 'color', 'y');
plot(target_set.slice([3,4], slice_at_vx_vy), 'color', 'k');
scatter(initial_state(1),initial_state(2),200,'ms','filled');
scatter([initial_state(1), optimal_mean_trajectory_ft(1,:)],...
        [initial_state(2), optimal_mean_trajectory_ft(2,:)],...
        30, 'bo', 'filled');
scatter([initial_state(1), optimal_mean_trajectory_cc_iter(1,:)],...
        [initial_state(2), optimal_mean_trajectory_cc_iter(2,:)],...
        30, 'co', 'filled');
scatter([initial_state(1), optimal_mean_trajectory_cc_pwl(1,:)],...
        [initial_state(2), optimal_mean_trajectory_cc_pwl(2,:)],...
        30, 'ro', 'filled');
legend_cell = {'Safe set','Target set','Initial state','Optimal mean trajectory (FT)', 'Optimal mean trajectory (CC-Iter)', 'Optimal mean trajectory (CC-PWL)'};
legend(legend_cell);

figure();
box on;
hold on;
plot(safe_set.slice([3,4], slice_at_vx_vy), 'color', 'y');
plot(target_set.slice([3,4], slice_at_vx_vy), 'color', 'k');
scatter(initial_state(1),initial_state(2),200,'ms','filled');
scatter([initial_state(1), optimal_mean_trajectory_ft(1,:)],...
        [initial_state(2), optimal_mean_trajectory_ft(2,:)],...
        30, 'bo', 'filled');
scatter([initial_state(1), optimal_mean_trajectory_cc_iter(1,:)],...
        [initial_state(2), optimal_mean_trajectory_cc_iter(2,:)],...
        30, 'co', 'filled');
scatter([initial_state(1), optimal_mean_trajectory_cc_pwl(1,:)],...
        [initial_state(2), optimal_mean_trajectory_cc_pwl(2,:)],...
        30, 'mo', 'filled');
legend_cell = {'Safe set','Target set','Initial state','Optimal mean trajectory (FT)', 'Optimal mean trajectory (CC-Iter)', 'Optimal mean trajectory (CC-PWL)'};
legend(legend_cell);

% Ft
green_legend_updated = 0;
red_legend_updated = 0;
traj_indices_ft = floor(n_mcarlo_sims*rand(1,n_sims_to_plot));
for realization_index = traj_indices_ft
    % Check if the trajectory satisfies the reach-avoid objective
    if mcarlo_result_ft(realization_index)
        % Assign green triangle as the marker
        markerString = 'g^-';
    else
        % Assign red asterisk as the marker
        markerString = 'r*-';
    end
    % Create [x(t_1) x(t_2)... x(t_N)]
    reshaped_X_vector = reshape(concat_state_realization_ft(:,realization_index), sys.state_dim,[]);
    % This realization is to be plotted
    h = plot([initial_state(1), reshaped_X_vector(1,:)], ...
             [initial_state(2), reshaped_X_vector(2,:)], ...
             markerString, 'MarkerSize',10);
    % Update the legends if the first else, disable
    if strcmp(markerString,'g^-')
        if green_legend_updated
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        else
            green_legend_updated = 1;
            legend_cell{end+1} = 'Good trajectory (FT)';
        end
    elseif strcmp(markerString,'r*-')
        if red_legend_updated
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        else
            red_legend_updated = 1;
            legend_cell{end+1} = 'Bad trajectory (FT)';
        end
    end
end
% CC (Iterative risk allocation)
green_legend_updated = 0;
red_legend_updated = 0;
traj_indices_cc_iter = floor(n_mcarlo_sims*rand(1,n_sims_to_plot));
for realization_index = traj_indices_cc_iter
    % Check if the trajectory satisfies the reach-avoid objective
    if mcarlo_result_cc_iter(realization_index)
        % Assign green triangle as the marker
        markerString = 'g^-';
    else
        % Assign red asterisk as the marker
        markerString = 'r*-';
    end
    % Create [x(t_1) x(t_2)... x(t_N)]
    reshaped_X_vector = reshape(concat_state_realization_cc_iter(:,realization_index), sys.state_dim,[]);
    % This realization is to be plotted
    h = plot([initial_state(1), reshaped_X_vector(1,:)], ...
             [initial_state(2), reshaped_X_vector(2,:)], ...
             markerString, 'MarkerSize',10);
    % Update the legends if the first else, disable
    if strcmp(markerString,'g^-')
        if green_legend_updated
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        else
            green_legend_updated = 1;
            legend_cell{end+1} = 'Good trajectory (Iter)';
        end
    elseif strcmp(markerString,'r*-')
        if red_legend_updated
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        else
            red_legend_updated = 1;
            legend_cell{end+1} = 'Bad trajectory (Iter)';
        end
    end
end
% CC (Piecewise linear --- PWL)
green_legend_updated = 0;
red_legend_updated = 0;
traj_indices_cc_pwl = floor(n_mcarlo_sims*rand(1,n_sims_to_plot));
for realization_index = traj_indices_cc_pwl
    % Check if the trajectory satisfies the reach-avoid objective
    if mcarlo_result_cc_iter(realization_index)
        % Assign green triangle as the marker
        markerString = 'g^-';
    else
        % Assign red asterisk as the marker
        markerString = 'r*-';
    end
    % Create [x(t_1) x(t_2)... x(t_N)]
    reshaped_X_vector = reshape(concat_state_realization_cc_pwl(:,realization_index), sys.state_dim,[]);
    % This realization is to be plotted
    h = plot([initial_state(1), reshaped_X_vector(1,:)], ...
             [initial_state(2), reshaped_X_vector(2,:)], ...
             markerString, 'MarkerSize',10);
    % Update the legends if the first else, disable
    if strcmp(markerString,'g^-')
        if green_legend_updated
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        else
            green_legend_updated = 1;
            legend_cell{end+1} = 'Good trajectory (PWL)';
        end
    elseif strcmp(markerString,'r*-')
        if red_legend_updated
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        else
            red_legend_updated = 1;
            legend_cell{end+1} = 'Bad trajectory (PWL)';
        end
    end
end
legend(legend_cell, 'Location','South');
title(sprintf('Plot with %d Monte-Carlo simulations from the initial state',...
                  n_sims_to_plot));
box on;
grid on;
% lb_stoch_reach_avoid_ft = 0;
% reach_avoid_probability_mcarlo_ft = 0;
fprintf('FT: %1.3f | CC (Iter): %1.3f | CC (PWL): %1.3f\n',...
    lb_stoch_reach_avoid_ft,...
    lb_stoch_reach_avoid_cc_iter,...
    lb_stoch_reach_avoid_cc_pwl); 
fprintf('MC (%1.0e particles): %1.3f, %1.3f, %1.3f\n',...
    n_mcarlo_sims,...
    sum(mcarlo_result_ft)/n_mcarlo_sims, ...
    sum(mcarlo_result_cc_iter)/n_mcarlo_sims, ...
    sum(mcarlo_result_cc_pwl)/n_mcarlo_sims);
fprintf('Elapsed time: %1.3f, %1.3f, %1.3f seconds\n',...
    elapsed_time_ft, elapsed_time_cc_iter, elapsed_time_cc_pwl);
