clear
clc
close all

ft_run = 0;
cc_iter_run = 0;
cc_pwl_run = 1;
cc_pwl_closed_run = 1;
cc_pwl_closed_cvx_run = 1; % For initial_state of -1,-1 and time horizon of 5, fails
                             % Works if time horizon shifted to 6 OR initial state is (-.75,-.75)

%% Monte-Carlo simulation parameters
n_mcarlo_sims = 1e5;
n_sims_to_plot = 5;
max_state_violation_prob = 0.01;
max_input_violation_prob = 0.01;
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
desired_accuracy_cc_closed = 5e-3;
PSoptions = psoptimset('Display','off');

%% Generate matrices for optimal mean trajectory generation
[H_matrix, mean_X_sans_input, ~] =...
       getHmatMeanCovForXSansInput(sys, initial_state, time_horizon);
%% 
% Next, using SReachTools, we will compute an <https://doi.org/10.1109/LCSYS.2017.2716364 
% optimal open-controller and the associated reach-avoid probability>. This function 
% takes about few minutes to run.

if ft_run
    timer_ft = tic;
    disp('Genz-algorithm+patternsearch technique');
    [lb_stoch_reach_avoid_ft, optimal_input_vector_ft] = ...
                             getLowerBoundStochReachAvoid(sys,...
                                                          initial_state,...
                                                          target_tube,...
                                                          'genzps',...
                                                          [],...
                                                          desired_accuracy,...
                                                          PSoptions);  
    elapsed_time_ft = toc(timer_ft);
    if lb_stoch_reach_avoid_ft > 0
        % This function returns the concatenated state vector stacked columnwise
        concat_state_realization_ft = generateMonteCarloSims(...
                                       n_mcarlo_sims,...
                                       sys,...
                                       initial_state,...
                                       time_horizon,...
                                       optimal_input_vector_ft);
        % Check if the location is within the target_set or not
        mcarlo_result_ft = target_tube.contains([repmat(initial_state,1,n_mcarlo_sims);
                                                 concat_state_realization_ft]);
        % Optimal mean trajectory generation                         
        optimal_mean_X_ft = mean_X_sans_input + H_matrix * optimal_input_vector_ft;
        optimal_mean_trajectory_ft=reshape(optimal_mean_X_ft,sys.state_dim,[]);                                              
    end
end

if cc_iter_run
    %% CC (Iterative risk allocation)     
    timer_cc_iter = tic;
    disp('Iterative risk allocation technique');
    [lb_stoch_reach_avoid_cc_iter, optimal_input_vector_cc_iter] =...
                             getLowerBoundStochReachAvoid(sys,...
                                                          initial_state,...
                                                          target_tube,...
                                                          'ccciter',...
                                                          desired_accuracy);  
    elapsed_time_cc_iter = toc(timer_cc_iter);
    if lb_stoch_reach_avoid_cc_iter>0
        % This function returns the concatenated state vector stacked columnwise
        concat_state_realization_cc_iter = generateMonteCarloSims(...
                                       n_mcarlo_sims,...
                                       sys,...
                                       initial_state,...
                                       time_horizon,...
                                       optimal_input_vector_cc_iter);
        % Check if the location is within the target_set or not
        mcarlo_result_cc_iter = target_tube.contains([repmat(initial_state,1,n_mcarlo_sims);
                                                      concat_state_realization_cc_iter]);
        % Optimal mean trajectory generation                         
        optimal_mean_X_cc_iter = mean_X_sans_input + H_matrix * optimal_input_vector_cc_iter;
        optimal_mean_trajectory_cc_iter=reshape(optimal_mean_X_cc_iter,sys.state_dim,[]);
    end
end

%% CC (Linear program approach)     
if cc_pwl_run
    timer_cc_pwl = tic;
    disp('Piecewise-linear single shot risk allocation technique');
    [lb_stoch_reach_avoid_cc_pwl, optimal_input_vector_cc_pwl] =...
                             getLowerBoundStochReachAvoid(sys,...
                                                          initial_state,...
                                                          target_tube,...
                                                          'cccpwl',...
                                                          desired_accuracy);  
    elapsed_time_cc_pwl = toc(timer_cc_pwl);
    if lb_stoch_reach_avoid_cc_pwl > 0
        % This function returns the concatenated state vector stacked columnwise
        concat_state_realization_cc_pwl = generateMonteCarloSims(...
                                       n_mcarlo_sims,...
                                       sys,...
                                       initial_state,...
                                       time_horizon,...
                                       optimal_input_vector_cc_pwl);
        % Check if the location is within the target_set or not
        mcarlo_result_cc_pwl = target_tube.contains([repmat(initial_state,1,n_mcarlo_sims);
                                                     concat_state_realization_cc_pwl]);
        % Optimal mean trajectory generation                         
        optimal_mean_X_cc_pwl = mean_X_sans_input + H_matrix * optimal_input_vector_cc_pwl;
        optimal_mean_trajectory_cc_pwl=reshape(optimal_mean_X_cc_pwl,sys.state_dim,[]);
    end
end

%% CC_pwl with closed loop (Linear program approach)     
if cc_pwl_closed_run
    timer_cc_pwl_closed = tic;
    disp('Piecewise-linear single shot (closed-loop-dc) risk allocation technique');
    [lb_stoch_reach_avoid_cc_pwl_closed, optimal_input_vector_cc_pwl_closed, optimal_input_gain, input_satisfaction_prob, risk_alloc] =...
                             getLowerBoundStochReachAvoid(sys,...
                                                          initial_state,...
                                                          target_tube,...
                                                          'cccpwl-closed',...
                                                          desired_accuracy_cc_closed,...
                                                          max_state_violation_prob,...
                                                          max_input_violation_prob);  
    elapsed_time_cc_pwl_closed = toc(timer_cc_pwl_closed);            
    if lb_stoch_reach_avoid_cc_pwl_closed > 0
        % This function returns the concatenated state vector stacked columnwise
        [concat_state_realization_cc_pwl_closed,concat_disturb_realization_cc_pwl_closed] = generateMonteCarloSims(...
                                       n_mcarlo_sims,...
                                       sys,...
                                       initial_state,...
                                       time_horizon,...
                                       optimal_input_vector_cc_pwl_closed,...
                                       optimal_input_gain);

        % Check if the location is within the target_set or not
        mcarlo_result_cc_pwl_closed = target_tube.contains([repmat(initial_state,1,n_mcarlo_sims);
                                                            concat_state_realization_cc_pwl_closed]);
        
        % Check if the input is within the tolerance
        [concat_input_space_A, concat_input_space_b] = sys.getConcatInputSpace(time_horizon);
        mcarlo_result_cc_pwl_closed_input = any(concat_input_space_A * (optimal_input_gain * concat_disturb_realization_cc_pwl_closed + optimal_input_vector_cc_pwl_closed)<=concat_input_space_b);
        
        % Optimal mean trajectory generation                         
        optimal_mean_X_cc_pwl_closed = mean_X_sans_input + H_matrix * optimal_input_vector_cc_pwl_closed;
        optimal_mean_trajectory_cc_pwl_closed=reshape(optimal_mean_X_cc_pwl_closed,sys.state_dim,[]);
    end
end

%% CC_pwl with closed loop (Linear program approach) with fixed slack
if cc_pwl_closed_cvx_run
    timer_cc_pwl_closed_cvx = tic;
    disp('Piecewise-linear single shot (closed-loop-cvx) with slack variables set to zero');
    [lb_stoch_reach_avoid_cc_pwl_closed_cvx, optimal_input_vector_cc_pwl_closed_cvx, optimal_input_gain_cvx, input_satisfaction_prob_cvx, risk_alloc_cvx] =...
                             getLowerBoundStochReachAvoid(sys,...
                                                          initial_state,...
                                                          target_tube,...
                                                          'cccpwl-closed-cvx',...
                                                          desired_accuracy_cc_closed,...
                                                          max_state_violation_prob,...
                                                          max_input_violation_prob);  
    elapsed_time_cc_pwl_closed_cvx = toc(timer_cc_pwl_closed_cvx);            
    if lb_stoch_reach_avoid_cc_pwl_closed_cvx > 0
        % This function returns the concatenated state vector stacked columnwise
        [concat_state_realization_cc_pwl_closed_cvx,concat_disturb_realization_cc_pwl_closed_cvx] = generateMonteCarloSims(...
                                       n_mcarlo_sims,...
                                       sys,...
                                       initial_state,...
                                       time_horizon,...
                                       optimal_input_vector_cc_pwl_closed_cvx,...
                                       optimal_input_gain_cvx);

        % Check if the location is within the target_set or not
        mcarlo_result_cc_pwl_closed_cvx = target_tube.contains([repmat(initial_state,1,n_mcarlo_sims);
                                                            concat_state_realization_cc_pwl_closed_cvx]);
        % Check if the input is within the tolerance
        [concat_input_space_A, concat_input_space_b] = sys.getConcatInputSpace(time_horizon);
        mcarlo_result_cc_pwl_closed_cvx_input = any(concat_input_space_A * (optimal_input_gain_cvx * concat_disturb_realization_cc_pwl_closed_cvx + optimal_input_vector_cc_pwl_closed_cvx)<=concat_input_space_b);
        
        % Optimal mean trajectory generation                         
        optimal_mean_X_cc_pwl_closed_cvx = mean_X_sans_input + H_matrix * optimal_input_vector_cc_pwl_closed_cvx;
        optimal_mean_trajectory_cc_pwl_closed_cvx =reshape(optimal_mean_X_cc_pwl_closed_cvx,sys.state_dim,[]);
    end
end

%% Plotting
disp('Plotting and Monte-Carlo simulation-based validation');
figure(1);
clf
box on;
hold on;
plot(safe_set.slice([3,4], slice_at_vx_vy), 'color', 'y');
plot(target_set.slice([3,4], slice_at_vx_vy), 'color', 'g');
scatter(initial_state(1),initial_state(2),200,'k^');
legend_cell = {'Safe set','Target set','Initial state'};
if exist('optimal_mean_trajectory_ft','var')
    scatter([initial_state(1), optimal_mean_trajectory_ft(1,:)],...
            [initial_state(2), optimal_mean_trajectory_ft(2,:)],...
            30, 'ro', 'filled');
    legend_cell{end+1} = 'Optimal mean trajectory (FT)';
end
if exist('optimal_mean_trajectory_cc_iter','var')
    scatter([initial_state(1), optimal_mean_trajectory_cc_iter(1,:)],...
        [initial_state(2), optimal_mean_trajectory_cc_iter(2,:)],...
        30, 'co', 'filled');
    legend_cell{end+1} = 'Optimal mean trajectory (CC-Iter)';
end
if exist('optimal_mean_trajectory_cc_pwl','var')
    scatter([initial_state(1), optimal_mean_trajectory_cc_pwl(1,:)],...
        [initial_state(2), optimal_mean_trajectory_cc_pwl(2,:)],...
        30, 'mo', 'filled');
    legend_cell{end+1} = 'Optimal mean trajectory (CC-PWL)';
end
if exist('optimal_mean_trajectory_cc_pwl_closed','var')
    scatter([initial_state(1), optimal_mean_trajectory_cc_pwl_closed(1,:)],...
        [initial_state(2), optimal_mean_trajectory_cc_pwl_closed(2,:)],...
        30, 'bo', 'filled');
    legend_cell{end+1} = 'Optimal mean trajectory (CC-PWL-Closed-DC)';
end
if exist('optimal_mean_trajectory_cc_pwl_closed_cvx','var')
    scatter([initial_state(1), optimal_mean_trajectory_cc_pwl_closed_cvx(1,:)],...
        [initial_state(2), optimal_mean_trajectory_cc_pwl_closed_cvx(2,:)],...
        30, 'ko', 'filled');
    legend_cell{end+1} = 'Optimal mean trajectory (CC-PWL-Closed-Cvx)';
end
legend(legend_cell);

%% Monte Carlo plot
figure(2);
clf
box on;
hold on;
plot(safe_set.slice([3,4], slice_at_vx_vy), 'color', 'y');
plot(target_set.slice([3,4], slice_at_vx_vy), 'color', 'g');
scatter(initial_state(1),initial_state(2),200,'k^');
legend_cell = {'Safe set','Target set','Initial state'};
if exist('optimal_mean_trajectory_ft','var')
    scatter([initial_state(1), optimal_mean_trajectory_ft(1,:)],...
            [initial_state(2), optimal_mean_trajectory_ft(2,:)],...
            30, 'ro', 'filled');
    legend_cell{end+1} = 'Optimal mean trajectory (FT)';
    if ~isnan(concat_state_realization_ft)
    %     [legend_cell] = plotMonteCarlo('(FT)', mcarlo_result_ft, concat_state_realization_ft, n_mcarlo_sims, n_sims_to_plot, sys.state_dim, initial_state, legend_cell);
        ellipsoidsFromMonteCarloSims(concat_state_realization_ft, sys.state_dim, [1,2], {'r'});
    end
else
    lb_stoch_reach_avoid_ft = NaN;
    mcarlo_result_ft = NaN;
    elapsed_time_ft = NaN;     
end
if exist('optimal_mean_trajectory_cc_iter','var')
    scatter([initial_state(1), optimal_mean_trajectory_cc_iter(1,:)],...
        [initial_state(2), optimal_mean_trajectory_cc_iter(2,:)],...
        30, 'co', 'filled');
    legend_cell{end+1} = 'Optimal mean trajectory (CC-Iter)';
    if ~isnan(concat_state_realization_cc_iter)
    %     [legend_cell] = plotMonteCarlo('(CC-Iter)', mcarlo_result_cc_iter, concat_state_realization_cc_iter, n_mcarlo_sims, n_sims_to_plot, sys.state_dim, initial_state, legend_cell);
        ellipsoidsFromMonteCarloSims(concat_state_realization_cc_iter, sys.state_dim, [1,2], {'c'});
    end
else
    lb_stoch_reach_avoid_cc_iter = NaN;
    mcarlo_result_cc_iter = NaN;
    elapsed_time_cc_iter = NaN;     
end
if exist('optimal_mean_trajectory_cc_pwl','var')
    scatter([initial_state(1), optimal_mean_trajectory_cc_pwl(1,:)],...
        [initial_state(2), optimal_mean_trajectory_cc_pwl(2,:)],...
        30, 'mo', 'filled');
    legend_cell{end+1} = 'Optimal mean trajectory (CC-PWL)';
    if ~isnan(concat_state_realization_cc_pwl)
%     [legend_cell] = plotMonteCarlo('(CC-PWL)', mcarlo_result_cc_pwl, concat_state_realization_cc_pwl, n_mcarlo_sims, n_sims_to_plot, sys.state_dim, initial_state, legend_cell);
        ellipsoidsFromMonteCarloSims(concat_state_realization_cc_pwl, sys.state_dim, [1,2], {'m'});
    end
else
    lb_stoch_reach_avoid_cc_pwl = NaN;
    mcarlo_result_cc_pwl = NaN;
    elapsed_time_cc_pwl = NaN;     
end
if exist('optimal_mean_trajectory_cc_pwl_closed','var')
        scatter([initial_state(1), optimal_mean_trajectory_cc_pwl_closed(1,:)],...
        [initial_state(2), optimal_mean_trajectory_cc_pwl_closed(2,:)],...
        30, 'bd', 'filled');
    legend_cell{end+1} = 'Optimal mean trajectory (CC-PWL-Closed-DC)';
    if ~isnan(concat_state_realization_cc_pwl_closed)
%         [legend_cell] = plotMonteCarlo('(CC-PWL-Closed-DC)', mcarlo_result_cc_pwl_closed, concat_state_realization_cc_pwl_closed, n_mcarlo_sims, n_sims_to_plot, sys.state_dim, initial_state, legend_cell);
        ellipsoidsFromMonteCarloSims(concat_state_realization_cc_pwl_closed, sys.state_dim, [1,2], {'b'});
    end
else
    lb_stoch_reach_avoid_cc_pwl_closed = NaN;
    mcarlo_result_cc_pwl_closed = NaN;
    elapsed_time_cc_pwl_closed = NaN;     
end
if exist('optimal_mean_trajectory_cc_pwl_closed_cvx','var')
    scatter([initial_state(1), optimal_mean_trajectory_cc_pwl_closed_cvx(1,:)],...
        [initial_state(2), optimal_mean_trajectory_cc_pwl_closed_cvx(2,:)],...
        30, 'ro', 'filled');
    legend_cell{end+1} = 'Optimal mean trajectory (CC-PWL-Closed-CVX)';
    if ~isnan(concat_state_realization_cc_pwl_closed_cvx)
    %     [legend_cell] = plotMonteCarlo('(CC-PWL-Closed-CVX)', mcarlo_result_cc_pwl_closed, concat_state_realization_cc_pwl_closed, n_mcarlo_sims, n_sims_to_plot, sys.state_dim, initial_state, legend_cell);
        ellipsoidsFromMonteCarloSims(concat_state_realization_cc_pwl_closed_cvx, sys.state_dim, [1,2], {'r'});
    end
else
    lb_stoch_reach_avoid_cc_pwl_closed_cvx = NaN;
    mcarlo_result_cc_pwl_closed_cvx = NaN;
    elapsed_time_cc_pwl_closed_cvx = NaN;     
end
legend(legend_cell, 'Location','South');
% title(sprintf('Plot with %d Monte-Carlo simulations from the initial state',...
%                   n_sims_to_plot));
title('Plot with ellipsoids fitting 100 randomly chosen Monte-Carlo simulations from the initial state');
box on;
grid on;
disp('Skipped items would show up as NaN');
fprintf('FT: %1.3f | CC (Iter): %1.3f | CC (PWL): %1.3f | CC (PWL-Closed-DC): %1.3f | CC (PWL-Closed-Cvx): %1.3f\n',...
    lb_stoch_reach_avoid_ft,...
    lb_stoch_reach_avoid_cc_iter,...
    lb_stoch_reach_avoid_cc_pwl,...
    lb_stoch_reach_avoid_cc_pwl_closed,...
    lb_stoch_reach_avoid_cc_pwl_closed_cvx); 
fprintf('MC (%1.0e particles): %1.3f, %1.3f, %1.3f, %1.3f, %1.3f\n',...
    n_mcarlo_sims,...
    sum(mcarlo_result_ft)/n_mcarlo_sims, ...
    sum(mcarlo_result_cc_iter)/n_mcarlo_sims, ...
    sum(mcarlo_result_cc_pwl)/n_mcarlo_sims,...
    sum(mcarlo_result_cc_pwl_closed)/n_mcarlo_sims,...
    sum(mcarlo_result_cc_pwl_closed_cvx)/n_mcarlo_sims);
fprintf('Elapsed time: %1.3f, %1.3f, %1.3f, %1.3f, %1.3f seconds\n',...
    elapsed_time_ft, elapsed_time_cc_iter, elapsed_time_cc_pwl, elapsed_time_cc_pwl_closed, elapsed_time_cc_pwl_closed_cvx);

function [legend_cell] = plotMonteCarlo(methodstr, mcarlo_result, concat_state_realization, n_mcarlo_sims, n_sims_to_plot, state_dim, initial_state, legend_cell)
    % CC (Piecewise linear --- PWL)
    green_legend_updated = 0;
    red_legend_updated = 0;
    traj_indices = floor(n_mcarlo_sims*rand(1,n_sims_to_plot));
    for realization_index = traj_indices
        % Check if the trajectory satisfies the reach-avoid objective
        if mcarlo_result(realization_index)
            % Assign green triangle as the marker
            markerString = 'g^-';
        else
            % Assign red asterisk as the marker
            markerString = 'r*-';
        end
        % Create [x(t_1) x(t_2)... x(t_N)]
        reshaped_X_vector = reshape(concat_state_realization(:,realization_index), state_dim,[]);
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
                legend_cell{end+1} = strcat('Good trajectory ', methodstr);
            end
        elseif strcmp(markerString,'r*-')
            if red_legend_updated
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            else
                red_legend_updated = 1;
                legend_cell{end+1} = strcat('Bad trajectory ', methodstr);
            end
        end
    end
end
