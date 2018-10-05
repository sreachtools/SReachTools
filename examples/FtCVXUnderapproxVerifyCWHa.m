clear
clc
close all

ft_run = 1;
cc_open_run = 1;
cc_affine_run = 1;

%% Monte-Carlo simulation parameters
n_mcarlo_sims = 1e5;
n_sims_to_plot = 5;
max_input_viol_prob = 0.01;
%% CWH system params
umax=0.1;
mean_disturbance = zeros(4,1);
covariance_disturbance = diag([1e-4, 1e-4, 5e-8, 5e-8]);
% Define the CWH (planar) dynamics of the deputy spacecraft relative to the 
% chief spacecraft as a LtiSystem object
sys = getCwhLtiSystem(4,...
                      Polyhedron('lb', -umax*ones(2,1),...
                                 'ub',  umax*ones(2,1)),...
                      StochasticDisturbance('Gaussian',...
                                            mean_disturbance,...
                                            covariance_disturbance));
time_horizon=5;          % Stay within a line of sight cone for 4 time steps and 
                         % reach the target at t=5% Safe Set --- LoS cone
%% Safe set definition --- LoS cone |x|<=y and y\in[0,ymax] and |vx|<=vxmax and 
%% |vy|<=vymax
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

%% Initial state definition
initial_state = [-1.15;         % Initial x relative position
                 -1.15;         % Initial y relative position
                 0;             % Initial x relative velocity
                 0];            % Initial y relative velocity
slice_at_vx_vy = initial_state(3:4);             
%% Parameters for MATLAB's Global Optimization Toolbox patternsearch
desired_accuracy = 1e-3;        % Decrease for a more accurate lower 
                                % bound at the cost of higher 
                                % computation time
desired_accuracy_cc_affine = 1e-2;
PSoptions = psoptimset('Display','off');

%% Generate matrices for optimal mean trajectory generation
[H_matrix, mean_X_sans_input, ~] =...
       getHmatMeanCovForXSansInput(sys, initial_state, time_horizon);

if ft_run
    timer_ft = tic;
    disp('Genz-algorithm+patternsearch technique');
    options.prob_str = 'term';
    options.method_str = 'genzps-open';
    options.desired_accuracy = desired_accuracy;
    options.PSoptions = psoptimset('display','off');
    [lb_stoch_reach_avoid_ft, optimal_input_vector_ft] = SReachPoint(...
        'term','genzps-open', sys, initial_state, target_tube, options);  
    elapsed_time_ft = toc(timer_ft);
    if lb_stoch_reach_avoid_ft > 0
        % This function returns the concatenated state vector stacked columnwise
        concat_state_realization_ft = generateMonteCarloSims(n_mcarlo_sims,...
            sys, initial_state, time_horizon, optimal_input_vector_ft);
        % Check if the location is within the target_set or not
        mcarlo_result_ft = target_tube.contains(...
            [repmat(initial_state,1,n_mcarlo_sims);
             concat_state_realization_ft]);
        % Optimal mean trajectory generation                         
        optimal_mean_X_ft = mean_X_sans_input +H_matrix*optimal_input_vector_ft;
        optimal_mean_trajectory_ft=reshape(optimal_mean_X_ft,sys.state_dim,[]);                                              
    end
end

%% CC (Linear program approach)     
if cc_open_run
    timer_cc_pwl = tic;
    disp('Piecewise-linear single shot risk allocation technique');
    options.prob_str = 'term';
    options.method_str = 'chance-open';
    options.desired_accuracy = desired_accuracy;
    [lb_stoch_reach_avoid_cc_pwl, optimal_input_vector_cc_pwl] = SReachPoint(...
        'term','chance-open', sys, initial_state, target_tube, options);  
    elapsed_time_cc_pwl = toc(timer_cc_pwl);
    if lb_stoch_reach_avoid_cc_pwl > 0
        % This function returns the concatenated state vector stacked columnwise
        concat_state_realization_cc_pwl = generateMonteCarloSims(...
            n_mcarlo_sims, sys, initial_state, time_horizon,...
            optimal_input_vector_cc_pwl);
        % Check if the location is within the target_set or not
        mcarlo_result_cc_pwl = target_tube.contains(...
            [repmat(initial_state,1,n_mcarlo_sims);
             concat_state_realization_cc_pwl]);
        % Optimal mean trajectory generation                         
        optimal_mean_X_cc_pwl = mean_X_sans_input +...
            H_matrix * optimal_input_vector_cc_pwl;
        optimal_mean_trajectory_cc_pwl=reshape(optimal_mean_X_cc_pwl,...
            sys.state_dim,[]);
    end
end

%% CC with affine controllers (Second order cone program approach)     
if cc_affine_run
    timer_cc_pwl_closed = tic;
    disp('Piecewise-linear single shot (closed-loop-dc) risk allocation technique');
    options.prob_str = 'term';
    options.method_str = 'chance-affine';
    options.desired_accuracy = desired_accuracy_cc_affine;
    options.max_input_viol_prob = max_input_viol_prob;
    options.bisect_lb = 0;
    options.bisect_ub = 1;
    options.tau_initial = 1;     % Weight for the DC constraint viol
    options.scaling_tau = 2;     % Multiplying factor for the weight
    options.tau_max = 1e5;       % Saturation for the weight to avoid numerical issues
    options.iter_max = 20;       % Max number of DC iterations
    options.dc_conv_tol = 1e-4;  % DC convergence threshold
    options.slack_tol = 1e-8;    % When is the slack variable == to the norm values?
    options.verbose = 1;

    [lb_stoch_reach_avoid_cc_pwl_closed, optimal_input_vector_cc_pwl_closed, optimal_input_gain, risk_alloc_state, risk_alloc_input] =...
         SReachPoint('term','chance-affine', sys, initial_state, target_tube,...
            options);  
    elapsed_time_cc_pwl_closed = toc(timer_cc_pwl_closed);            
    if lb_stoch_reach_avoid_cc_pwl_closed > 0
        % This function returns the concatenated state vector stacked columnwise
        [concat_state_realization_cc_pwl_closed,...
            concat_disturb_realization_cc_pwl_closed] =...
                generateMonteCarloSims(n_mcarlo_sims, sys, initial_state,...
                    time_horizon,optimal_input_vector_cc_pwl_closed,...
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
    legend_cell{end+1} = 'Optimal mean trajectory (genzps-open)';
end
if exist('optimal_mean_trajectory_cc_pwl','var')
    scatter([initial_state(1), optimal_mean_trajectory_cc_pwl(1,:)],...
        [initial_state(2), optimal_mean_trajectory_cc_pwl(2,:)],...
        30, 'mo', 'filled');
    legend_cell{end+1} = 'Optimal mean trajectory (chance-open)';
end
if exist('optimal_mean_trajectory_cc_pwl_closed','var')
    scatter([initial_state(1), optimal_mean_trajectory_cc_pwl_closed(1,:)],...
        [initial_state(2), optimal_mean_trajectory_cc_pwl_closed(2,:)],...
        30, 'bo', 'filled');
    legend_cell{end+1} = 'Optimal mean trajectory (chance-affine)';
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
    legend_cell{end+1} = 'Optimal mean trajectory (genzps-open)';
    if ~isnan(concat_state_realization_ft)
%         [legend_cell] = plotMonteCarlo('(genzps-open)', mcarlo_result_ft,...
%             concat_state_realization_ft, n_mcarlo_sims, n_sims_to_plot,...
%             sys.state_dim, initial_state, legend_cell);
        ellipsoidsFromMonteCarloSims(concat_state_realization_ft, sys.state_dim, [1,2], {'r'});
    end
else
    lb_stoch_reach_avoid_ft = NaN;
    mcarlo_result_ft = NaN;
    elapsed_time_ft = NaN;     
end
if exist('optimal_mean_trajectory_cc_pwl','var')
    scatter([initial_state(1), optimal_mean_trajectory_cc_pwl(1,:)],...
        [initial_state(2), optimal_mean_trajectory_cc_pwl(2,:)],...
        30, 'mo', 'filled');
    legend_cell{end+1} = 'Optimal mean trajectory (chance-open)';
    if ~isnan(concat_state_realization_cc_pwl)
%         [legend_cell] = plotMonteCarlo('(chance-open)', ...
%             mcarlo_result_cc_pwl, concat_state_realization_cc_pwl,...
%             n_mcarlo_sims, n_sims_to_plot, sys.state_dim, initial_state,...
%             legend_cell);
        ellipsoidsFromMonteCarloSims(concat_state_realization_cc_pwl,...
            sys.state_dim, [1,2], {'m'});
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
    legend_cell{end+1} = 'Optimal mean trajectory (chance-affine)';
    if ~isnan(concat_state_realization_cc_pwl_closed)
%         [legend_cell] = plotMonteCarlo('(chance-affine)',...
%             mcarlo_result_cc_pwl_closed,...
%             concat_state_realization_cc_pwl_closed, n_mcarlo_sims,...
%             n_sims_to_plot, sys.state_dim, initial_state, legend_cell);
        ellipsoidsFromMonteCarloSims(concat_state_realization_cc_pwl_closed,...
            sys.state_dim, [1,2], {'b'});
    end
else
    lb_stoch_reach_avoid_cc_pwl_closed = NaN;
    mcarlo_result_cc_pwl_closed = NaN;
    elapsed_time_cc_pwl_closed = NaN;     
end
legend(legend_cell, 'Location','South');
% title(sprintf('Plot with %d Monte-Carlo sims', n_sims_to_plot));
title('Plot with ellipsoid fit for 100 randomly chosen Monte-Carlo sims');
box on;
grid on;
disp('Skipped items would show up as NaN');
fprintf('FT: %1.3f | CC (Open): %1.3f | CC (Affine): %1.3f\n',...
    lb_stoch_reach_avoid_ft,...
    lb_stoch_reach_avoid_cc_pwl,...
    lb_stoch_reach_avoid_cc_pwl_closed); 
fprintf('MC (%1.0e particles): %1.3f, %1.3f, %1.3f\n',...
    n_mcarlo_sims,...
    sum(mcarlo_result_ft)/n_mcarlo_sims, ...
    sum(mcarlo_result_cc_pwl)/n_mcarlo_sims,...
    sum(mcarlo_result_cc_pwl_closed)/n_mcarlo_sims);
fprintf('Elapsed time: %1.3f, %1.3f, %1.3f seconds\n',...
    elapsed_time_ft, elapsed_time_cc_pwl, elapsed_time_cc_pwl_closed);

%%
function [legend_cell] = plotMonteCarlo(method_str, mcarlo_result,...
    concat_state_realization, n_mcarlo_sims, n_sims_to_plot, state_dim,...
    initial_state, legend_cell)
% Plots a selection of Monte-Carlo simulations on top of the plot

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
        reshaped_X_vector = reshape(...
            concat_state_realization(:,realization_index), state_dim,[]);
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
                legend_cell{end+1} = strcat('Good trajectory ', method_str);
            end
        elseif strcmp(markerString,'r*-')
            if red_legend_updated
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            else
                red_legend_updated = 1;
                legend_cell{end+1} = strcat('Bad trajectory ', method_str);
            end
        end
    end
end
