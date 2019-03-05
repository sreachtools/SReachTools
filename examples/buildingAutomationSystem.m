%% Building automation system verification
% See ARCH 2019
%
% Comparison of (chance-constraint+open-loop)-based underapproximation with
% Lagrangian (set-theoretic)-based underapproximation to construct the
% stochastic reach set at prob_thresh = 0.8
%
% In the interest of time, Genzps+patternsearch-based underapproximation has
% been disabled. See should_we_run_genzps

clear;clc;close all;srtinit;

%% Problem setup
% System matrices
state_matrix = ...
   [ 0.9096   , -0.0106  , -0.0019   ,        0;
    -0.0097   , 0.9152  ,       0   , -0.0015;
    -0.0133   , 0        ,  0.9625  ,      0;
    0         , -0.0133  ,       0   , 0.9625];
input_matrix = [0.0779;0.0779;0;0];
dist_matrix = eye(4);

% Input space
input_space = Polyhedron('lb',19,'ub',20);

% Disturbance definition
dist_mu = [0.0622;0.0561;1.8094;1.8094];
dist_sigma = diag([0.2,0.2,0.1,0.1]);
dist_rv = RandomVector.gaussian(dist_mu, dist_sigma);

% System definition
sys = LtiSystem('StateMatrix',state_matrix, 'InputMatrix',input_matrix, ...
    'DisturbanceMatrix',dist_matrix, 'InputSpace', input_space, ...
    'Disturbance', dist_rv);

% Time steps
time_horizon = 6;

% Safety specification --- Constraints are present only the first two states
safe_set = Polyhedron('lb', [-19.5;-19.5;-Inf;-Inf], 'ub', [20.5;20.5;Inf;Inf]);
safety_tube = Tube('viability', safe_set, time_horizon);

% Stochastic viability threshold: Compute the set of initial states which
% have the probability of safety above this threshold
prob_thresh = 0.8;
% Slice of the stochastic viability set of interest
x3_init = 0;
x4_init = 0;

% How many directions to explore for chance-open and genzps
cc_n_dir_vecs = 16;
genzps_n_dir_vecs = 16;
lag_over_n_dir_vecs = 16;
should_we_run_genzps = 0; % DISABLED: It takes time for (slightly) better result

%% Construction of the chance-open-based underapproximation
disp('>>> Chance-constraint-based underapproximation');
timerVal = tic;
% Directions to explore
theta_vec = linspace(0, 2*pi, cc_n_dir_vecs + 1);
theta_vec = theta_vec(1:end-1);
set_of_dir_vecs = [cos(theta_vec);sin(theta_vec);zeros(2,cc_n_dir_vecs)];
% Slice of the stochastic viability set of interest
init_safe_set_affine = Polyhedron('He',[0,0,1,0,x3_init;0,0,0,1,x4_init]);
% SReachSet options preparation
cco_options = SReachSetOptions('term', 'chance-open', ...
    'init_safe_set_affine', init_safe_set_affine, 'set_of_dir_vecs', ...
    set_of_dir_vecs, 'verbose', 1);
cco_stoch_viab_set = SReachSet('term','chance-open',sys, prob_thresh, ...
    safety_tube, cco_options);
elapsed_time_cc = toc(timerVal);
% For plotting, construct the slice
cco_stoch_viab_set_2D =  cco_stoch_viab_set.slice([3,4], [x3_init;x4_init]);

%% Construction of the lagrangian-based underapproximation
fprintf('\n\n\n >>> Lagrangian-based underapproximation\n');
timerVal = tic;
lag_under_options = SReachSetOptions('term', 'lag-under', ...
    'bound_set_method', 'ellipsoid', 'compute_style', 'vfmethod', 'verbose', 1);
lag_under_stoch_viab_set = SReachSet('term','lag-under',sys, prob_thresh, ...
    safety_tube, lag_under_options);
elapsed_time_lag_under = toc(timerVal);
% For plotting, construct the slice
lag_under_stoch_viab_set_2D = lag_under_stoch_viab_set.slice([3,4], ...
    [x3_init;x4_init]);

% %% Construction of the lagrangian-based overapproximation
% fprintf('\n\n\n >>> Lagrangian-based overapproximation\n');
% timerVal = tic;
% lag_over_options = SReachSetOptions('term', 'lag-over', ...
%     'bound_set_method', 'ellipsoid', 'compute_style', 'support', ...
%     'verbose', 1, 'system', sys, 'n_vertices', lag_over_n_dir_vecs);
% lag_over_stoch_viab_set = SReachSet('term','lag-over',sys, prob_thresh, ...
%     safety_tube, lag_over_options);
% elapsed_time_lag_over = toc(timerVal);
% % For plotting, construct the slice
% lag_over_stoch_viab_set_2D = lag_over_stoch_viab_set.slice([3,4], ...
%     [x3_init;x4_init]);

%% Construction of the Genz+patternsearch-based underapproximation
if should_we_run_genzps
    fprintf('\n\n\n >>> Genz+patternsearch-based underapproximation\n');
    timerVal = tic;
    % Directions to explore
    theta_vec = linspace(0, 2*pi, genzps_n_dir_vecs + 1);
    theta_vec = theta_vec(1:end-1);
    set_of_dir_vecs = [cos(theta_vec);sin(theta_vec);
                       zeros(2,genzps_n_dir_vecs)];
    % Slice of the stochastic viability set of interest
    init_safe_set_affine = Polyhedron('He',[0,0,1,0,x3_init;0,0,0,1,x4_init]);
    % SReachSet options preparation
    genzps_options = SReachSetOptions('term', 'genzps-open', ...
        'init_safe_set_affine', init_safe_set_affine, 'set_of_dir_vecs', ...
        set_of_dir_vecs,'verbose',1);
    genzps_stoch_viab_set = SReachSet('term','genzps-open',sys, prob_thresh, ...
        safety_tube, genzps_options);
    elapsed_time_genzps = toc(timerVal);
    % For plotting, construct the slice
    genzps_stoch_viab_set_2D =  genzps_stoch_viab_set.slice([3,4], ...
        [x3_init; x4_init]);
else
    fprintf('\n\n\n>>> Skipping genzps-based computation.\n');
end

%% Plot the figures
figure(1);
clf
plot(Polyhedron('V',safe_set.V(:,1:2)),'color','y');
hold on;
if should_we_run_genzps
    plot(genzps_stoch_viab_set_2D, 'color','g');
end
plot(cco_stoch_viab_set_2D, 'color','m');
plot(lag_under_stoch_viab_set_2D, 'color','b');
if should_we_run_genzps
    leg=legend('Safe set', 'Stochastic viability set (genzps-open)', ...
        'Stochastic viability set (chance-open)', ...
        'Stochastic viability set (lag-under)');
else
    leg=legend('Safe set', 'Stochastic viability set (chance-open)', ...
        'Stochastic viability set (lag-under)');
end
set(leg,'Location','EastOutside');
title(sprintf('Safety analysis for $x_3[0]=$%1.2f, $x_4[0]=$%1.2f', x3_init, ...
    x4_init), 'interpreter','latex');
box on;
grid on;
axis equal;
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
% In general, increase the fontsize
set(gca,'FontSize',20);
% If code ocean, save the results
% saveas(gcf, '../results/BAS_StochasticViabilitySet.png');

%% Disp
if should_we_run_genzps
    fprintf(['\n\n\n>>> Compute times\n Chance constraint-based ',...
        'underapprox.: %1.2f\n Lagrangian-based underapprox.: %1.2f\n ', ...
        'Genzps-based computation: %1.2f\n'], ...
        elapsed_time_cc, elapsed_time_lag_under, elapsed_time_genzps);
else
    fprintf(['\n\n\n>>> Compute times\n Chance constraint-based ',...
        'underapprox.: %1.2f\n Lagrangian-based underapprox.: %1.2f\n'], ...
        elapsed_time_cc, elapsed_time_lag_under);
end

