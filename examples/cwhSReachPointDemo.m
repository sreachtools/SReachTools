%% Demonstration of controller synthesis via SReachPoint: Spacecraft rendezvous
% This example will demonstrate the use of SReachTools in controller synthesis
% for stochastic continuous-state discrete-time linear time-invariant (LTI) 
% systems.
% 
% Specifically, we will discuss how we can use SReachPoint to synthesize 
% open-loop controllers and affine-disturbance feedback controllers for
% that maximize the probability of safety while respecting the system dynamics 
% and control bounds. We demonstrate:
% 
% * |chance-open|: Chance-constrained approach that uses risk allocation and 
%    piecewise-affine approximations to formulate a linear program to
%    synthesize an open-loop controller
% * |genzps-open|: Fourier transforms that uses
%    <http://www.math.wsu.edu/faculty/genz/software/matlab/qsimvnv.m Genz's 
%    algorithm> to formulate a nonlinear log-concave optimization problem to be
%    solved using MATLAB's patternsearch to synthesize an open-loop controller
% * |particle-open|: Particle control approach filter that formulates a 
%    mixed-integer linear program to synthesize an open-loop controller
% * |chance-affine|: Chance-constrained approach that uses risk allocation and 
%    piecewise-affine approximations to formulate a difference-of-convex program 
%    to synthesize a closed-loop (affine disturbance feedback) controller.
%    The controller synthesis is done by solving a series of second-order
%    cone programs.
%
% Our approaches are grid-free and recursion-free resulting in highly scalable
% solutions, especially for Gaussian-perturbed LTI systems. 
%
% This script is part of the SReachTools toolbox, which is licensed under GPL v3
% or (at your option) any later version. A copy of this license is given in 
% <https://github.com/unm-hscl/SReachTools/blob/master/LICENSE 
% https://github.com/unm-hscl/SReachTools/blob/master/LICENSE>.

% Prescript running
close all;
% clc;
clearvars;
srtinit
%% Problem formulation: Spacecraft motion via CWH dynamics
% We consider both the spacecrafts, referred to as the deputy spacecraft and 
% the chief spacecraft, to be in the same circular orbit. In this example, we 
% will consider the problem of optimal controller synthesis for the deputy such 
% that it can rendezvous with the chief while staying within the line-of-sight
% cone with maximium likelihood.
%%
% <<cwh_sketch.png>>
%% Dynamics model for the deputy relative to the chief spacecraft
% The relative planar dynamics of the deputy with respect to the chief are described 
% by the <https://doi.org/10.1109/CDC.2013.6760626 Clohessy-Wiltshire-Hill (CWH) 
% equations,> 
% 
% $$\ddot{x} - 3 \omega x - 2 \omega \dot{y} = \frac{F_{x}}{m_{d}}$$
% 
% $$            \ddot{y} + 2 \omega \dot{x} = \frac{F_{y}}{m_{d}}$$ 
% 
% where the position of the deputy relative to the chief is $x,y \in
% \mathbf{R}$, $\omega = \sqrt{\frac{\mu}{R_{0}^{3}}}$ is the orbital frequency,
% $\mu$ is the gravitational constant, and $R_{0}$ is the orbital radius of the
% chief spacecraft.  We define the state as $\overline{x} = {[x\ y\ \dot{x}\
% \dot{y}]}^\top \in \mathbf{R}^{4}$ which is the position and velocity of the
% deputy relative to the chief along $\mathrm{x}$- and $\mathrm{y}$- axes, and
% the input as $\overline{u} = {[F_{x}\ F_{y}]}^\top \in
% \mathcal{U}\subset\mathbf{R}^{2}$. 
% 
% We will discretize the CWH dynamics in time, via zero-order hold, to obtain
% the discrete-time linear time-invariant system and add a Gaussian disturbance
% to account for the modeling uncertainties and the disturbance forces,
% 
% $$\overline{x}_{k+1} = A \overline{x}_{k} + B \overline{u}_{k} +
% \overline{w}_{k}$$
% 
% with $\overline{w}_{k} \in \mathbf{R}^{4}$ as an IID Gaussian zero-mean 
% random process with a known covariance matrix $\Sigma_{\overline{w}}$. 
 
%% System definition
umax = 0.1;
mean_disturbance = zeros(4,1);
covariance_disturbance = diag([1e-4, 1e-4, 5e-8, 5e-8]);
% Define the CWH (planar) dynamics of the deputy spacecraft relative to the
% chief spacecraft as a LtiSystem object
sys = getCwhLtiSystem(4, Polyhedron('lb', -umax*ones(2,1), ...
                                        'ub',  umax*ones(2,1)), ...
       RandomVector('Gaussian', mean_disturbance,covariance_disturbance));

%% Target tube construction --- reach-avoid specification
time_horizon = 5;   % Stay within a line of sight cone for 4 time steps and 
                    % reach the target at t=5% Safe Set --- LoS cone
% Safe set definition --- LoS cone |x|<=y and y\in[0,ymax] and |vx|<=vxmax and 
% |vy|<=vymax
ymax = 2;
vxmax = 0.5;
vymax = 0.5;
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
% Target set --- Box [-0.1,0.1]x[-0.1,0]x[-0.01,0.01]x[-0.01,0.01]
target_set = Polyhedron('lb', [-0.1; -0.1; -0.01; -0.01], ...
                        'ub', [0.1; 0; 0.01; 0.01]);
target_tube = Tube('reach-avoid',safe_set, target_set, time_horizon);


%% Specifying initial states and which options to run
cc_open_run = 1;
ft_run = 1;
pa_open_run = 1;
cc_affine_run = 1;
% Initial state definition
initial_state = [-1.15;         % Initial x relative position
                 -1.15;         % Initial y relative position
                 0;             % Initial x relative velocity
                 0];            % Initial y relative velocity
slice_at_vx_vy = initial_state(3:4); 
init_state_ccc_open = initial_state;
init_state_genzps_open = initial_state;
init_state_particle_open = initial_state;
init_state_ccc_affine = initial_state;


%% Preparation for Monte-Carlo simulations of the optimal controllers
% Monte-Carlo simulation parameters
n_mcarlo_sims = 1e5;
n_sims_to_plot = 5;      % Required only if plot_traj_instead_of_ellipses = 1
% Generate matrices for optimal mean trajectory generation
% Get H and mean_X_sans_input
[~, H, ~] = getConcatMats(sys, time_horizon);
sysnoi = LtvSystem('StateMatrix',sys.state_mat,'DisturbanceMatrix', ...
    sys.dist_mat,'Disturbance',sys.dist);
[mean_X_sans_input, ~] = SReachFwd('concat-stoch', sysnoi, initial_state, ...
    time_horizon);

%% SReachPoint: |chance-open|
% This implements the chance-constrained approach to compute a globally opt 
% open-loop controller. This approach uses risk allocation and piecewise-affine 
% overapproximation of the inverse normal cumulative density function to
% formulate a linear program for this purpose. Naturally, this is one of
% the fastest ways to compute an open-loop controller and an
% underapproximative probabilistic guarantee of safety. However, due to the use
% of Boole's inequality for risk allocation, it provides a conservative estimate
% of safety using the open-loop controller.
if cc_open_run
    timer_cc_pwl = tic;
    [lb_stoch_reach_avoid_cc_pwl, optimal_input_vector_cc_pwl] = SReachPoint(...
        'term','chance-open', sys, initial_state, target_tube);  
    elapsed_time_cc_pwl = toc(timer_cc_pwl);
    if lb_stoch_reach_avoid_cc_pwl > 0
        % This function returns the concatenated state vector stacked columnwise
        concat_state_realization_cc_pwl = generateMonteCarloSims(...
            n_mcarlo_sims, sys, initial_state, time_horizon, ...
            optimal_input_vector_cc_pwl);
        % Check if the location is within the target_set or not
        mcarlo_result_cc_pwl = target_tube.contains(...
            concat_state_realization_cc_pwl);
        % Optimal mean trajectory generation                         
        optimal_mean_X_cc_pwl = mean_X_sans_input + ...
            H * optimal_input_vector_cc_pwl;
        optimal_mean_trajectory_cc_pwl = reshape(optimal_mean_X_cc_pwl, ...
            sys.state_dim,[]);
    end
end

%% SReachPoint: |genzps-open|
% This implements the Fourier transform-based approach to compute a globally
% opt open-loop controller. This approach uses
% <http://www.math.wsu.edu/faculty/genz/software/matlab/qsimvnv.m Genz's
% algorithm> to compute the probability of safety and optimizes the joint chance
% constraint involved in maximizing this probability. To handle the noisy
% behaviour of the Genz's algorithm, we rely on MATLAB's |patternsearch| for the
% nonlinear optimization. The global optity of the open-loop controller is
% guaranteed by the log-concavity of the problem. Internally, we use the
% |chance-open| to initialize the nonlinear solver. Hence, this approach will
% return an open-loop controller with safety at least as good as |chance-open|.
if ft_run
    timer_ft = tic;
    [lb_stoch_reach_avoid_ft, optimal_input_vector_ft] = SReachPoint(...
        'term','genzps-open', sys, initial_state, target_tube);  
    elapsed_time_ft = toc(timer_ft);
    if lb_stoch_reach_avoid_ft > 0
        % This function returns the concatenated state vector stacked columnwise
        concat_state_realization_ft = generateMonteCarloSims(n_mcarlo_sims, ...
            sys, initial_state, time_horizon, optimal_input_vector_ft);
        % Check if the location is within the target_set or not
        mcarlo_result_ft = target_tube.contains(concat_state_realization_ft);
        % Optimal mean trajectory generation                         
        optimal_mean_X_ft = mean_X_sans_input + H * optimal_input_vector_ft;
        optimal_mean_trajectory_ft = reshape(optimal_mean_X_ft,sys.state_dim,[]);                                              
    end
end

%% SReachPoint: |chance-affine|
% This implements the chance-constrained approach to compute a locally opt 
% affine disturbance feedback controller. This approach uses risk allocation and 
% piecewise-affine overapproximation of the inverse normal cumulative density 
% function to formulate a difference-of-convex program. We utilize penalty 
% convex-concave procedure to solve this program to a local optimum. Due to
% its construction of the affine feedback controller, this approach
% typically permits the construction of the highest underapproximative
% probability guarantee. Since affine disturbance feedback controllers can
% not satisfy hard control bounds, we relax the control bounds to be
% probabilistically violated with at most a probability of 0.01
if cc_affine_run
    timer_cc_affine = tic;
    options = SReachPointOptions('term','chance-affine', ...
        'max_input_viol_prob', 1e-2, 'verbose', 1);
    [lb_stoch_reach_avoid_cc_affine, optimal_input_vector_cc_affine, ...
        optimal_input_gain, risk_alloc_state, risk_alloc_input] =...
         SReachPoint('term','chance-affine', sys, initial_state, target_tube, ...
            options);  
    elapsed_time_cc_affine = toc(timer_cc_affine);            
    if lb_stoch_reach_avoid_cc_affine > 0
        % This function returns the concatenated state vector stacked columnwise
        [concat_state_realization_cc_affine, ...
            concat_disturb_realization_cc_affine] =...
                generateMonteCarloSims(n_mcarlo_sims, sys, initial_state, ...
                    time_horizon,optimal_input_vector_cc_affine, ...
                    optimal_input_gain);

        % Check if the location is within the target_set or not
        mcarlo_result_cc_affine = target_tube.contains(...
            concat_state_realization_cc_affine);
        
        % Check if the input is within the tolerance
        [concat_input_space_A, concat_input_space_b] =...
            sys.getConcatInputSpace(time_horizon);
        mcarlo_result_cc_affine_input = any(concat_input_space_A * ...
            (optimal_input_gain * concat_disturb_realization_cc_affine + ...
                optimal_input_vector_cc_affine)<=concat_input_space_b);
        
        % Optimal mean trajectory generation                         
        optimal_mean_X_cc_affine = mean_X_sans_input + ...
            H * optimal_input_vector_cc_affine;
        optimal_mean_trajectory_cc_affine = reshape(...
            optimal_mean_X_cc_affine,sys.state_dim,[]);
    end
end

%% Particle filter approach
if pa_open_run
    timer_pa = tic;
    [lb_stoch_reach_avoid_pa, optimal_input_vector_pa] = SReachPoint(...
        'term','particle-open', sys, initial_state, target_tube);  
    elapsed_time_pa = toc(timer_pa);
    if lb_stoch_reach_avoid_pa > 0
        % This function returns the concatenated state vector stacked columnwise
        concat_state_realization_pa = generateMonteCarloSims(...
            n_mcarlo_sims, sys, initial_state, time_horizon, ...
            optimal_input_vector_pa);
        % Check if the location is within the target_set or not
        mcarlo_result_pa = target_tube.contains(concat_state_realization_pa);
        % Optimal mean trajectory generation                         
        optimal_mean_X_pa = mean_X_sans_input + ...
            H * optimal_input_vector_pa;
        optimal_mean_trajectory_pa = reshape(optimal_mean_X_pa, ...
            sys.state_dim,[]);
    end
end

%% Plotting and Monte-Carlo simulation-based validation
figure(1);
clf
box on;
hold on;
plot(safe_set.slice([3,4], slice_at_vx_vy), 'color', 'y');
plot(target_set.slice([3,4], slice_at_vx_vy), 'color', 'g');
scatter(initial_state(1),initial_state(2),200,'k^');
legend_cell = {'Safe set','Target set','Initial state'};
if exist('optimal_mean_trajectory_ft','var')
    scatter([initial_state(1), optimal_mean_trajectory_ft(1,:)], ...
            [initial_state(2), optimal_mean_trajectory_ft(2,:)], ...
            30, 'ro', 'filled');
    legend_cell{end+1} = 'Optimal mean trajectory (genzps-open)';
end
if exist('optimal_mean_trajectory_cc_pwl','var')
    scatter([initial_state(1), optimal_mean_trajectory_cc_pwl(1,:)], ...
        [initial_state(2), optimal_mean_trajectory_cc_pwl(2,:)], ...
        30, 'mo', 'filled');
    legend_cell{end+1} = 'Optimal mean trajectory (chance-open)';
end
if exist('optimal_mean_trajectory_pa','var')
    scatter([initial_state(1), optimal_mean_trajectory_pa(1,:)], ...
        [initial_state(2), optimal_mean_trajectory_pa(2,:)], ...
        30, 'ko', 'filled');
    legend_cell{end+1} = 'Optimal mean trajectory (particle-open)';
end
if exist('optimal_mean_trajectory_cc_affine','var')
    scatter([initial_state(1), optimal_mean_trajectory_cc_affine(1,:)], ...
        [initial_state(2), optimal_mean_trajectory_cc_affine(2,:)], ...
        30, 'bo', 'filled');
    legend_cell{end+1} = 'Optimal mean trajectory (chance-affine)';
end
legend(legend_cell, 'Location','South');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');

figure(2);
clf
box on;
hold on;
plot(safe_set.slice([3,4], slice_at_vx_vy), 'color', 'y');
plot(target_set.slice([3,4], slice_at_vx_vy), 'color', 'g');
scatter(initial_state(1),initial_state(2),200,'k^');
legend_cell = {'Safe set','Target set','Initial state'};
if exist('optimal_mean_trajectory_ft','var')
    scatter([initial_state(1), optimal_mean_trajectory_ft(1,:)], ...
            [initial_state(2), optimal_mean_trajectory_ft(2,:)], ...
            30, 'ro', 'filled');
    legend_cell{end+1} = 'Optimal mean trajectory (genzps-open)';
    if ~isnan(concat_state_realization_ft)
        if plot_traj_instead_of_ellipses == 1
            [legend_cell] = plotMonteCarlo('(genzps-open)', mcarlo_result_ft, ...
                concat_state_realization_ft, n_mcarlo_sims, n_sims_to_plot, ...
                sys.state_dim, initial_state, legend_cell);
        else
            ellipsoidsFromMonteCarloSims(concat_state_realization_ft, ...
                sys.state_dim, [1,2], {'r'});
        end
    end
else
    lb_stoch_reach_avoid_ft = NaN;
    mcarlo_result_ft = NaN;
    elapsed_time_ft = NaN;     
end
if exist('optimal_mean_trajectory_cc_pwl','var')
    scatter([initial_state(1), optimal_mean_trajectory_cc_pwl(1,:)], ...
        [initial_state(2), optimal_mean_trajectory_cc_pwl(2,:)], ...
        30, 'mo', 'filled');
    legend_cell{end+1} = 'Optimal mean trajectory (chance-open)';
    if ~isnan(concat_state_realization_cc_pwl)
        if plot_traj_instead_of_ellipses == 1
            [legend_cell] = plotMonteCarlo('(chance-open)', ...
                mcarlo_result_cc_pwl, concat_state_realization_cc_pwl, ...
                n_mcarlo_sims, n_sims_to_plot, sys.state_dim, initial_state, ...
                legend_cell);
        else
            ellipsoidsFromMonteCarloSims(concat_state_realization_cc_pwl, ...
                sys.state_dim, [1,2], {'m'});
        end
    end
else
    lb_stoch_reach_avoid_cc_pwl = NaN;
    mcarlo_result_cc_pwl = NaN;
    elapsed_time_cc_pwl = NaN;     
end
if exist('optimal_mean_trajectory_pa','var')
    scatter([initial_state(1), optimal_mean_trajectory_pa(1,:)], ...
        [initial_state(2), optimal_mean_trajectory_pa(2,:)], ...
        30, 'ko', 'filled');
    legend_cell{end+1} = 'Optimal mean trajectory (particle-open)';
    if ~isnan(concat_state_realization_pa)
        if plot_traj_instead_of_ellipses == 1
            [legend_cell] = plotMonteCarlo('(particle-open)', ...
                mcarlo_result_pa, concat_state_realization_pa, ...
                n_mcarlo_sims, n_sims_to_plot, sys.state_dim, initial_state, ...
                legend_cell);
        else
            ellipsoidsFromMonteCarloSims(concat_state_realization_pa, ...
                sys.state_dim, [1,2], {'k'});
        end
    end
else
    lb_stoch_reach_avoid_pa = NaN;
    mcarlo_result_pa = NaN;
    elapsed_time_pa = NaN;     
end

if exist('optimal_mean_trajectory_cc_affine','var')
        scatter([initial_state(1), optimal_mean_trajectory_cc_affine(1,:)], ...
        [initial_state(2), optimal_mean_trajectory_cc_affine(2,:)], ...
        30, 'bd', 'filled');
    legend_cell{end+1} = 'Optimal mean trajectory (chance-affine)';
    if ~isnan(concat_state_realization_cc_affine)
        if plot_traj_instead_of_ellipses==1
            [legend_cell] = plotMonteCarlo('(chance-affine)', ...
                mcarlo_result_cc_affine, ...
                concat_state_realization_cc_affine, n_mcarlo_sims, ...
                n_sims_to_plot, sys.state_dim, initial_state, legend_cell);
        else
            ellipsoidsFromMonteCarloSims(concat_state_realization_cc_affine, ...
                sys.state_dim, [1,2], {'b'});
        end
    end
else
    lb_stoch_reach_avoid_cc_affine = NaN;
    mcarlo_result_cc_affine = NaN;
    elapsed_time_cc_affine = NaN;     
end
legend(legend_cell, 'Location','South');
if plot_traj_instead_of_ellipses==1
    title(sprintf('Plot with %d Monte-Carlo sims', n_sims_to_plot));
else
    title('Plot with ellipsoid fit for 100 randomly chosen Monte-Carlo sims');
end
box on;
grid on;
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
%% Reporting the results
if any(isnan([lb_stoch_reach_avoid_cc_pwl, lb_stoch_reach_avoid_ft, ...
        lb_stoch_reach_avoid_cc_affine, lb_stoch_reach_avoid_pa]))
    disp('Skipped items would show up as NaN');
end
fprintf(['FT: %1.3f | CC (Open): %1.3f | Scenario (Open): %1.3f | ', ...
    'CC (Affine): %1.3f\n'], ...
    lb_stoch_reach_avoid_ft, ...
    lb_stoch_reach_avoid_cc_pwl, ...
    lb_stoch_reach_avoid_pa, ...
    lb_stoch_reach_avoid_cc_affine); 
fprintf('MC (%1.0e particles): %1.3f, %1.3f, %1.3f, %1.3f\n', ...
    n_mcarlo_sims, ...
    sum(mcarlo_result_ft)/n_mcarlo_sims, ...
    sum(mcarlo_result_cc_pwl)/n_mcarlo_sims, ...
    sum(mcarlo_result_pa)/n_mcarlo_sims, ...
    sum(mcarlo_result_cc_affine)/n_mcarlo_sims);
fprintf('Elapsed time: %1.3f, %1.3f, %1.3f, %1.3f seconds\n', ...
    elapsed_time_ft, elapsed_time_cc_pwl, elapsed_time_pa, ...
    elapsed_time_cc_affine);

%% Helper functions
function [legend_cell] = plotMonteCarlo(method_str, mcarlo_result, ...
    concat_state_realization, n_mcarlo_sims, n_sims_to_plot, state_dim, ...
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
