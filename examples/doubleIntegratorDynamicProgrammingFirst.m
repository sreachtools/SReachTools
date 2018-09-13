%% Double Integrator Reach-Avoid Via Dynamic Programming
% This example demonstrates how to use the SReachTools toolbox to solve a terminal-hitting 
% time reach-avoid problem using <https://doi.org/10.1016/j.automatica.2010.08.006 
% dynamic programming>.
% 
% In this example, we analyze the following problems via dynamic programming 
% for a stochastic system with known dynamics:
% 
% # *stochastic viability problem: *Compute a controller to stay within a safe 
% set with maximum likelihood
% # *the terminal-hitting time stochastic reach-avoid problem: *Compute a controller 
% that maximizes the probability of reaching a target set at a time horizon, |N|, 
% while maintaining the system in a set of safe states
% # *stochastic reachability of a moving target tube: *Compute a controller 
% that maximizes the probability of staying within a target tube
% 
% SReachTools has a dynamic programing implementation that can analyze systems 
% upto three dimensions. For efficient implementation, we require the input set 
% to be an axis-aligned hypercuboid, and define the grid the smallest hypercuboid 
% containing all the target sets.
%% Notes about this Live Script:
% # *MATLAB dependencies*: This Live Script uses MATLAB's  <https://www.mathworks.com/products/statistics.html 
% Statistics and Machine Learning Toolbox>. 
% # *External dependencies*: This Live Script uses Multi-Parameteric Toolbox 
% (<http://people.ee.ethz.ch/~mpt/3/ MPT>). 
% # Make sure that |srtinit| is run before running this script.
% 
% This Live Script is part of the SReachTools toolbox. License for the use 
% of this function is given in <https://github.com/unm-hscl/SReachTools/blob/master/LICENSE 
% https://github.com/unm-hscl/SReachTools/blob/master/LICENSE>.
%% Problem setup
%% Double Integrator
% In this example we use a discretized double integrator dynamics given by:
% 
% $$  x_{k+1} = \left[ \begin{array}{cc}    1 & T \\    0 & 1  \end{array}\right] 
% x_{k} + \left[\begin{array}{c}    \frac{T^{2}}{2} \\    T  \end{array}\right] 
% u_{k} + w_{k}$$
% 
% where $<math xmlns="http://www.w3.org/1998/Math/MathML" display="inline"><mrow><mi 
% mathvariant="italic">T</mi></mrow></math>$ is the discretization time-step, 
% and $w_{k}$ is the stochastic disturbance.
%% Setup the system
%%
% discretization parameter
T = 0.25;

% define the system
sys = LtiSystem('StateMatrix', [1, T; 0, 1], ...
    'InputMatrix', [T^2/2; T], ...
    'InputSpace', Polyhedron('lb', -0.1, 'ub', 0.1), ...
    'DisturbanceMatrix', eye(2), ...
    'Disturbance', StochasticDisturbance('Gaussian', zeros(2,1), 0.001*eye(2)));
%% Setup the dynamic programming and visualization parameters

dyn_prog_xinc = 0.05;
dyn_prog_uinc = 0.1;
reach_set_thresholds = [0.2 0.5 0.9];
legend_str={'Target tube at t=0', 'Safety Probability $\geq 0.2$', 'Safety Probability $\geq 0.5$', 'Safety Probability $\geq 0.9$'};
%% Case 1: Terminal hitting-time stochatsic reach-avoid problem (Constant target sets uptil time horizon - 1 and potentially a different target set at time horizon)
%% Setup the target and safe sets
%%
safe_set = Polyhedron('lb', [-1, -1], 'ub', [1, 1]);
target_set = Polyhedron('V',[-.2,0;1,1;1,-1]);%Polyhedron('lb', -[0.5, 0.5], 'ub', [0.5, 0.5]);
axis_vec1 = [-1 1 -1 1];
%% Setup the target tube
% Target tube is a generalization of the reach problem. The reach avoid target-tube 
% is created by setting the first $<math xmlns="http://www.w3.org/1998/Math/MathML" 
% display="inline"><mrow><mi mathvariant="italic">N</mi><mo>?</mo><mn>1</mn></mrow></math>$ 
% sets in the tube as the |safe_set| and the final set as the |target_set|.

% time horizon
N = 15;
% in target tube for the viability problem is equivalent to a tube of repeating
% safe sets
safety_tube1 = TargetTube('viability', safe_set, N);

% Plotting of target tube
figure()
hold on    
for time_indx=0:N
    target_tube_at_time_indx = Polyhedron('H',[safety_tube1(time_indx+1).A,zeros(size(safety_tube1(time_indx+1).A,1),1), safety_tube1(time_indx+1).b], 'He',[0 0 1 time_indx]);
    plot(target_tube_at_time_indx, 'alpha',0.25);
end
axis([axis_vec1 0 N])
box on;
grid on;
xlabel('x');
ylabel('y');
zlabel('time');
title('Target tube');
%% Dynamic programming recursion via gridding
%%
tic;
[prob_x1, cell_of_xvec_x1] = getDynProgSolForTargetTubeFirst(sys, dyn_prog_xinc, dyn_prog_uinc, safety_tube1,target_set);
toc
%% Visualization of the value function at t=0 (safety probability)
%%
figure();
x1vec = cell_of_xvec_x1{1};
x2vec = cell_of_xvec_x1{2};
surf(x1vec,x2vec,reshape(prob_x1,length(x2vec),length(x1vec)));
axis([axis_vec1 0 1])
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
zlabel('Safety probability')
box on
view(45, 45)
%% Visualization of the safe initial sets --- Superlevel sets of safety probability

figure();
poly_array1 = getDynProgLevelSets2D(cell_of_xvec_x1, prob_x1, reach_set_thresholds, safety_tube1);
hold on;
plot([safety_tube1(1), poly_array1])
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
box on
axis(axis_vec1)
axis equal
legend(legend_str,'interpreter','latex');


%% Case 2: Hyperplane target set
N=5;
safe_set2 = Polyhedron('lb', [-1, -1], 'ub', [1, 1]);
% -x -y <- 1 => x + y > 1
target_hyperplane = Polyhedron('H',[-1 -1 -1]);
axis_vec2 = [-1 1 -1 1];
safety_tube2 = TargetTube('viability', safe_set2, N);
tic
[prob_x2, cell_of_xvec_x2, grid_x] = getDynProgSolForTargetTubeFirst(sys, dyn_prog_xinc/2, dyn_prog_uinc, safety_tube2, target_hyperplane);
toc
% getLowerBoundStochReachAvoidFirst(sys,...
%                                 initial_state,...
%                                 safety_tube,...
%                                 target_hyperplane,...
%                                 varargin)

%% Plot it
figure();
x1vec = cell_of_xvec_x2{1};
x2vec = cell_of_xvec_x2{2};
surf(x1vec,x2vec,reshape(prob_x2,length(x2vec),length(x1vec)));
axis([axis_vec2 0 1 ])
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
zlabel('Safety probability')
box on
view(45, 45)

%% Chance constrained underapproximation
l_xmax = 0.15;
xvec_for_plot = l_xmax+dyn_prog_xinc/2:dyn_prog_xinc/2:1;
dyn_soln = prob_x2(min(grid_x')>=l_xmax);
useful_grid_x = grid_x(min(grid_x')>=l_xmax,:);
lb_stoch_reach = zeros(length(useful_grid_x),1);
tic;
for ix = 1:length(useful_grid_x)
    initial_state = useful_grid_x(ix,:)';
    [lb_stoch_reach(ix)]=getLowerBoundStochReachAvoidFirst(sys,...
                                    initial_state,...
                                    safety_tube2,...
                                    target_hyperplane,...
                                    1e-3);
    disp([ix/length(useful_grid_x) sum(initial_state)-1 dyn_soln(ix)-lb_stoch_reach(ix)]);    
end
toc
%% Plot
figure()
surf(xvec_for_plot,xvec_for_plot,reshape(lb_stoch_reach,length(xvec_for_plot),[]))
hold on;
surf(xvec_for_plot,xvec_for_plot,reshape(dyn_soln,length(xvec_for_plot),[]))
