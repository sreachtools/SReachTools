% Prescript running: Initializing srtinit, if it already hasn't been initialized
close all;clearvars;srtinit;srtinit --version;

time_horizon = 50;
time_const = 1/2*time_horizon;
init_heading = pi/10;
sampling_time = 0.1;
box_halflength = 4;
omega = pi/time_horizon/sampling_time;
turning_rate = omega*ones(time_horizon,1);
dist_cov = 0.0001;
prob_thresh = 0.99;
no_of_direction_vectors_ccc = 16;
v_nominal = 10;
umax = v_nominal/3*2;
n_mcarlo_sims = 1e3;

%% LTV system definition
[sys, heading_vec] = getDubinsCarLtv('add-dist', ...
turning_rate, ...
init_heading, ...
sampling_time, ...
Polyhedron('lb',0,'ub',umax), ...
eye(2), ...
RandomVector('Gaussian',zeros(2,1), dist_cov * eye(2)));

target_tube_cell = cell(time_horizon + 1,1);

%% Target tube definition
figure(100);clf;hold on
angle_at_the_center = (heading_vec) - pi/2;
center_box = zeros(2, time_horizon + 1);        
for itt = 0:time_horizon
    center_box(:, itt+1) = v_nominal * [cos(angle_at_the_center(itt+1))-cos(angle_at_the_center(1));
                                        sin(angle_at_the_center(itt+1))-sin(angle_at_the_center(1))];
    target_tube_cell{itt+1} = Polyhedron('lb',center_box(:, itt+1) - box_halflength * exp(- itt/time_const), 'ub', center_box(:, itt+1) + box_halflength*exp(- itt/time_const));
    plot(target_tube_cell{itt+1},'alpha',0.5,'color','y');
end
axis equal
axis([-8    10   -5   21]);
box on;
grid on;

target_tube = Tube(target_tube_cell{:});

%% Set of direction vectors
theta_vector_ccc = linspace(0, 2*pi, no_of_direction_vectors_ccc+1);
theta_vector_ccc = theta_vector_ccc(1:end-1);
set_of_direction_vectors_ccc = [cos(theta_vector_ccc); 
                                sin(theta_vector_ccc)];
% %% Set computation              
% timer_polytope_ccc = tic;
% opts = SReachSetOptions('term', 'chance-open', 'pwa_accuracy', 1e-3, ...
%         'set_of_dir_vecs', set_of_direction_vectors_ccc, ...
%         'init_safe_set_affine',Polyhedron(),'verbose',0);
% [ccc_polytope, extra_info] = SReachSet('term','chance-open', sys, ...
%       prob_thresh, target_tube, opts);
% elapsed_time_polytope_ccc = toc(timer_polytope_ccc);
% fprintf('Time taken for computing the polytope (CCC): %1.3f s\n', elapsed_time_polytope_ccc);

%% Lagrangian under
timer_lagunder = tic;
theta_polytope_vec = linspace(0,2*pi,10)';
lagunder_options = SReachSetOptions('term', 'lag-under', ...
    'bound_set_method', 'polytope', 'template_polytope', ...
    Polyhedron('V',[cos(theta_polytope_vec),sin(theta_polytope_vec)]), ...
    'compute_style', 'vfmethod', 'vf_enum_method', 'lrs', 'verbose', 2);
%     'bound_set_method', 'ellipsoid', 'compute_style', 'vfmethod', ...
%     'vf_enum','lrs', 'verbose', 2);

[polytope_lagunder, extra_info_under] = SReachSet('term', 'lag-under', ...
    sys, prob_thresh, target_tube, lagunder_options);
elapsed_time_lagunder = toc(timer_lagunder);

% Compute a far-away safe initial
cvx_begin quiet
    variable initial_state(sys.state_dim, 1)
    minimize ([1 1]*initial_state)
    subject to
        polytope_lagunder.A*initial_state <= polytope_lagunder.b;
        target_tube(1).A*initial_state <= target_tube(1).b;
cvx_end
switch cvx_status
    case 'Solved'
        fprintf('Testing initial state: ');
        disp(initial_state');
        
        % Create a controller based on the underapproximation
        srlcontrol = SReachLagController(sys, ... 
            extra_info_under.bounded_dist_set, ...
            extra_info_under.stoch_reach_tube);
        % Generate Monte-Carlo simulations using the srlcontrol and
        % generateMonteCarloSims
        timer_mcarlo = tic;
        [X,U,W] = generateMonteCarloSims(n_mcarlo_sims, sys, ...
            initial_state, time_horizon, srlcontrol, [], ...
            lagunder_options.verbose);
        elapsed_time_mcarlo = toc(timer_mcarlo);
        avg_time_mc = elapsed_time_mcarlo / n_mcarlo_sims;

        % % Plot the convex hull of the spread of the points
        polytopesFromMonteCarloSims(X, 4, [1,2], {'color','k','alpha',0});
        a = gca;            
        for tindx = 1:time_horizon-1
            a.Children(tindx).Annotation.LegendInformation.IconDisplayStyle='off';
        end
        a.Children(1).Annotation.LegendInformation.IconDisplayStyle='on';
        a.Children(1).DisplayName = 'Trajectory spread at various time steps';           
        
        % Plot the initial state
        scatter(initial_state(1), initial_state(2), 200, 'ko', 'filled', ...
            'DisplayName','Initial state');
    otherwise        
end

init_state_lag = initial_state;
optimal_mean_trajectory_lag = reshape(sum(X, 2) / size(X, 2), 2, []);


%% Plot the set
figure(101);
clf;
hold on;
plot(target_tube(1));
% plot(underapproximate_stochastic_reach_avoid_polytope_ccc,'color','m');
plot(ccc_polytope,'color','b');
axis equal
axis (1.2*[-box_halflength box_halflength -box_halflength box_halflength]);
box on;
legend('Target set at t=0','Stochastic reach set','Location','SouthEast');
set(gca,'FontSize',20);
