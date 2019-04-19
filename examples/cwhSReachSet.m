%% Verification of satellite rendezvous problem via SReachSet 
% This example will demonstrate the use of |SReachTools| for controller
% synthesis in a terminal hitting-time stochastic reach-avoid problem. We
% consider a continuous-state discrete-time linear time-invariant (LTI) system.
% This example script is part of the |SReachTools| toolbox, which is licensed
% under GPL v3 or (at your option) any later version. A copy of this license is
% given in <https://sreachtools.github.io/license/
% https://sreachtools.github.io/license/>.
% 
% Specifically, we will discuss how SReachTools can use Fourier transforms
% (<http://www.math.wsu.edu/faculty/genz/software/matlab/qsimvnv.m Genz's
% algorithm> and MATLAB's patternsearch), convex chance constraints, and
% Lagrangian methods to construct underapproximative stochastic reach sets.
% Our approaches are grid-free and recursion-free, resulting in highly scalable
% solutions. 
%
% All computations were performed using MATLAB on an Ubuntu OS running on a
% laptop with Intel i7 CPU with 2.1GHz clock rate and 8 GB RAM. For sake of
% clarity, all commands were asked to be verbose (via `SReachSetOptions`). In
% practice, this can be turned off.

% Prescript running: Initializing srtinit, if it already hasn't been initialized
close all;clearvars;srtinit;

%% Problem formulation: Spacecraft motion via CWH dynamics
% We consider both the spacecrafts, referred to as the deputy spacecraft and 
% the chief spacecraft, to be in the same circular orbit. In this example, we 
% will consider the problem of verification for the spacecraft rendezvous
% problem, i.e., identify all the initial states from which the deputy can can
% rendezvous with the chief while staying within the line-of-sight cone with
% a likelihood above a user-specified threshold.
%%
% <<cwh_sketch.png>>
%% Dynamics model for the deputy relative to the chief spacecraft 
% The relative planar dynamics of the deputy with respect to the chief are
% described by the <https://doi.org/10.1109/CDC.2013.6760626
% Clohessy-Wiltshire-Hill (CWH) equations,> 
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
 
umax = 0.1;
mean_disturbance = zeros(4,1);
covariance_disturbance = diag([1e-4, 1e-4, 5e-8, 5e-8]);
% Define the CWH (planar) dynamics of the deputy spacecraft relative to the
% chief spacecraft as a LtiSystem object
sys = getCwhLtiSystem(4, Polyhedron('lb', -umax*ones(2,1), ...
                                    'ub',  umax*ones(2,1)), ...
       RandomVector('Gaussian', mean_disturbance,covariance_disturbance));

%% Target tube construction --- reach-avoid specification
% We define the target tube to be a collection of time-varying boxes
% $\{\mathcal{T}_k\}_{k=0}^N$ where $N$ is the time horizon.
%
% In this problem, we define $\mathcal{T}_k$ to be line-of-sight cone 
% originating from origin (location of the chief spacecraft) for
% $k\in\{0,1,\ldots,N-1\}$ and the terminal target set $\mathcal{T}_N$ as a
% box around the origin. This special sequence of target sets allows us to
% impose a reach-avoid specification of safety.

time_horizon = 5;        % Stay within a line of sight cone for 4 time steps and 
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


%% Preparation for set computation

% Specify restrictions on the initial velocity
slice_at_vx_vy = zeros(2,1);      
init_safe_set_affine = Polyhedron('He',[zeros(2,2) eye(2,2) slice_at_vx_vy]);
prob_thresh = 0.8;
% Direction vectors for Fourier transform approach
n_dir_vecs = 15;
theta_vec = linspace(0, 2*pi, n_dir_vecs);
set_of_dir_vecs_ft = [cos(theta_vec);
                      sin(theta_vec);
                      zeros(2,n_dir_vecs)];
% Direction vectors for chance-constrained approach
n_dir_vecs = 40;
theta_vec = linspace(0, 2*pi, n_dir_vecs);
set_of_dir_vecs_cc_open = [cos(theta_vec);
                           sin(theta_vec);
                           zeros(2,n_dir_vecs)];
% THIS SCRIPT ONLY --- Flags that enable running specific components
ft_run = 1;
cc_open_run = 1;
lagunder_run = 1;
lagover_run = 0;

%% Verification via |SReachSet|
% These approaches utilizes the convexity and compactness guarantees of the 
% stochastic reach set to the problem of stochastic reachability of a target
% tube as discussed in
%
% # A. Vinod and M. Oishi, "Scalable underapproximative verification of stochastic
% LTI systems using convexity and compactness," In Proc. Hybrid Syst.: Comput. &
% Ctrl., pages 1--10, 2018. HSCC 2018
% # A. Vinod and M. Oishi, "Stochastic reachability of a target tube: Theory and
% computation," IEEE Transactions in Automatic Control, 2018 (submitted)
% https://arxiv.org/pdf/1810.05217.pdf.
%
% The three different approaches explored in this example are
%
% # Chance-constrained open-loop-based verification (Linear program approach)
% # Genz's algorithm+MATLAB's patternsearch+open-loop-based verification
% # Lagrangian-based underapproximation
%
% While the first two methods use ray-shooting and |SReachPoint| to compute
% a polytopic underapproximation, the third approach utilizes
% Lagrangian-based underapproximation as described in
%
% * J. D. Gleason, A. P. Vinod, M. M. K. Oishi, "Underapproximation of
%   Reach-Avoid Sets for Discrete-Time Stochastic Systems via Lagrangian
%   Methods," in Proceedings of the IEEE Conference on Decision and Control,
%   2017.

% Convex chance constrained-based underapproximation
if cc_open_run
    cc_options = SReachSetOptions('term', 'chance-open', ...
        'set_of_dir_vecs', set_of_dir_vecs_cc_open, ...
        'init_safe_set_affine', init_safe_set_affine, ...
        'verbose', 1);
    timer_cc_open = tic;
    [polytope_cc_open, extra_info] = SReachSet('term','chance-open', sys, ...
        prob_thresh, target_tube, cc_options);  
    elapsed_time_cc_open = toc(timer_cc_open);
end
%%

% Genz's algorithm + Patternsearch-based underapproximation
if ft_run
    ft_options = SReachSetOptions('term', 'genzps-open', ...
        'set_of_dir_vecs', set_of_dir_vecs_ft, ...
        'init_safe_set_affine', init_safe_set_affine, 'verbose', 1, ...
        'desired_accuracy', 5e-2);
    timer_ft = tic;
    polytope_ft = SReachSet('term','genzps-open', sys, prob_thresh, ...
        target_tube, ft_options);  
    elapsed_time_ft = toc(timer_ft);
end
%%

% Lagrangian-based underapproximation
if lagunder_run
    n_dim = sys.state_dim + sys.input_dim;
    lagunder_options = SReachSetOptions('term', 'lag-under',...
        'bound_set_method', 'ellipsoid', 'compute_style','support',...
        'system', sys, 'vf_enum_method', 'lrs', 'verbose', 2,...
        'n_vertices', 2^n_dim * 6 + 2*n_dim);
    
    timer_lagunder = tic;
    [polytope_lagunder, extra_info_under] = SReachSet('term', 'lag-under', ...
        sys, prob_thresh, target_tube, lagunder_options);
    elapsed_time_lagunder = toc(timer_lagunder);
end
%%

% Lagrangian-based overapproximation
if lagover_run
    n_dim = sys.state_dim;
    lagover_options = SReachSetOptions('term', 'lag-over', ...
        'bound_set_method', 'ellipsoid', 'compute_style','support', ...
        'system', sys, 'vf_enum_method', 'lrs', 'verbose', 1, ...
        'n_vertices', 2^n_dim * 6 + 2*n_dim);
    
    timer_lagover = tic;
    polytope_lagover = SReachSet('term', 'lag-over', sys,  prob_thresh, ...
        target_tube, lagover_options);
    elapsed_time_lagover = toc(timer_lagover);
end

%% Plotting the sets
figure(101);
clf
box on;
hold on;
plot(safe_set.slice([3,4], slice_at_vx_vy), 'color', 'y');
plot(target_set.slice([3,4], slice_at_vx_vy), 'color', 'g');
legend_cell = {'Safe set','Target set'};
if cc_open_run && ~isEmptySet(polytope_cc_open)
    % Since we specify the init_safe_set (polytope is already intersected by 
    % the constant vx and vy at t=0), we can simply throw out the
    % states 3 and 4 to get the 2D set
    plot(Polyhedron('V',polytope_cc_open.V(:,1:2)), 'color','k','alpha',0.5);
    legend_cell{end+1} = 'Underapprox. polytope (chance-open)';
else
    polytope_cc_open = Polyhedron();
    elapsed_time_cc_open = NaN;
end
if ft_run && ~isEmptySet(polytope_ft)
    % Since we specify the init_safe_set (polytope is already intersected by 
    % the constant vx and vy at t=0), we can simply throw out the
    % states 3 and 4 to get the 2D set
    plot(Polyhedron('V',polytope_ft.V(:,1:2)), 'color','b','alpha',0.5);
    legend_cell{end+1} = 'Underapprox. polytope (genzps-open)';
else
    polytope_ft = Polyhedron();
    elapsed_time_ft = NaN;
end
if lagunder_run && ~isEmptySet(polytope_lagunder)
    % Since we specify the init_safe_set (polytope is already intersected by 
    % the constant vx and vy at t=0), we can simply throw out the
    % states 3 and 4 to get the 2D set
    plot(polytope_lagunder.slice([3,4], slice_at_vx_vy), 'color','r','alpha',1);
    legend_cell{end+1} = 'Underapprox. polytope (lag-under)';
else
    polytope_lagunder = Polyhedron();
    elapsed_time_lagunder = NaN;
end
if lagover_run && ~isEmptySet(polytope_lagover)
    plot(polytope_lagover.slice([3,4], slice_at_vx_vy), 'color','m','alpha',1);
    legend_cell{end+1} = 'Overapprox. polytope (lag-over)';
else
    polytope_lagover = Polyhedron();
    elapsed_time_lagover = NaN;
end

%% Monte-Carlo based validation of the optimal controller

% Monte-Carlo simulation parameters
n_mcarlo_sims = 1e5;
n_sims_to_plot = 5;
direction_index_to_plot = 30;
if ~isEmptySet(polytope_cc_open)
    init_state = extra_info(2).vertices_underapprox_polytope(:,direction_index_to_plot);
    input_vec = extra_info(2).opt_input_vec_at_vertices(:,direction_index_to_plot);
    opt_reach_avoid = extra_info(2).opt_reach_prob_i(direction_index_to_plot);

    concat_state_realization = generateMonteCarloSims(...
            n_mcarlo_sims, ...
            sys, ...
            init_state, ...
            time_horizon, ...
            input_vec);        
    
    % Check if the location is within the target_set or not
    mcarlo_result = target_tube.contains(concat_state_realization);
    [legend_cell] = plotMonteCarlo(' (chance-open)', mcarlo_result, ...
        concat_state_realization, n_mcarlo_sims, n_sims_to_plot, ...
        sys.state_dim, init_state, legend_cell);
    fprintf('Expected probability: %1.3f, Simulated probability: %1.3f\n',...
        opt_reach_avoid, sum(mcarlo_result)/n_mcarlo_sims);

end
legend(legend_cell, 'Location','South');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
if any(isnan([elapsed_time_ft, elapsed_time_cc_open, elapsed_time_lagunder]))
    disp('Skipped items would show up as NaN');
end
fprintf(['Elapsed time: (genzps-open) %1.3f | (chance-open) %1.3f',...
    ' | (lag-under) %1.3f | (lag-over) %1.3f seconds\n'], ...
    elapsed_time_ft, elapsed_time_cc_open, elapsed_time_lagunder,...
    elapsed_time_lagover);

%%
%% A helper function to plot a subset of the Monte-Carlo simulations
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
