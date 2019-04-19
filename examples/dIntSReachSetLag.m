%% Lagrangian Approximations for the Stochastic Reachability of a Target Tube
% This example will demonstrate the use of |SReachTools| for Lagrangian-based 
% verification of stochastic continuous-state discrete-time linear
% time-invariant (LTI) systems.  This example script is part of the
% |SReachTools| toolbox, which is licensed under GPL v3 or (at your option) any
% later version. A copy of this license is given in
% <https://sreachtools.github.io/license/
% https://sreachtools.github.io/license/>.
% 
%% Lagrangian Methods
% Lagrangian methods perform computations with sets using operations like unions, 
% intersection, Minkowski addition/differences, etc. This computation using set 
% operations can be used to approximate (either over or under) the stochastic 
% reachability of a target tube problem. We will demonstrate that this approach, 
% while being be approximative, can outperform the current state-of-the-art
% <https://doi.org/10.1016/j.automatica.2010.08.006 dynamic programming>
% solution in terms of computation time.
% 
% Advantages:
% 
% * No gridding, which partially evades the curse of dimensionality
% * Provides verification for closed-loop feedback strategies
% * Synthesis of a closed-loop feedback strategy
% 
% Disadvantages:
% 
% * Using Polyhedral representation, must solve the vertex-facet enumeration 
% problem, limiting computations to ~4 dimensional systems
% 
% The theory for this approach can be found in
% 
% * J. D. Gleason, A. P. Vinod, M. M. K. Oishi, "Underapproximation of
%   Reach-Avoid Sets for Discrete-Time Stochastic Systems via Lagrangian
%   Methods," in Proceedings of the IEEE Conference on Decision and Control,
%   2017. 
%
% Further, we explore multiple implementations of the Lagrangian-based
% verification, where the vertex-facet enumeration is mitigated either via
% recursion-free support method in overapproximation or vertex-complexity
% preserving support vector method in underapproximation.
%
% All computations were performed using MATLAB on an Ubuntu OS running on a
% laptop with Intel i7 CPU with 2.1GHz clock rate and 8 GB RAM. For sake of
% clarity, all commands were asked to be verbose (via `SReachSetOptions`). In
% practice, this can be turned off.

% Prescript running: Initializing srtinit, if it already hasn't been initialized
close all;clearvars;srtinit;

%% Problem Definition
% In this example we will look at the viability problem for a double integrator. 
% The system dynamics are:
% 
% $$  x_{k+1} = \left[ \begin{array}{cc}    1 & T \\   0 & 1  \end{array}\right] 
% x_{k} + \left[\begin{array}{c}    \frac{T^{2}}{2} \\    T  \end{array}\right] 
% u_{k} + w_{k}$$
% 
% where $x_{k+1} \in \mathbf{R}^{2}$ is the state, $u_{k} \in \mathbf{R}$ 
% is the input, and $w_{k} \in \mathbf{R}^{2}$ is the disturbance. The following 
% code defines this system with $w_{k}$ as an i.i.d. Gaussian disturbance with 
% mean $[0, 0]^{\top}$ and variance $\mathrm{diag}(0.01, 0.01)$.
%%
% example parameters
T = 0.25;

% define the system
sys = getChainOfIntegLtiSystem(2, ...
    T, ...
    Polyhedron('lb', -0.1, 'ub', 0.1), ...
     RandomVector('Gaussian', zeros(2,1), 0.001*eye(2)));

%% Viability problem as a stochastic reachability of a target tube problem
% We examine the viability problem in which we are interested in staying in 
% a set of safe states. In this example the safe set is $\{x \in \mathbf{R}^{2}: 
% |x_{i}| < 1, i = 1, 2\}$. The stochastic reachability of a target tube problem 
% posed as a viability problem by constructing a target tube in which all sets 
% in the tube are the safe set.
%%
time_horizon = 5;

% safe set definition
safe_set = Polyhedron('lb', [-1, -1], 'ub', [1, 1]);
% target tube definition
target_tube = Tube('viability', safe_set, time_horizon);
% probability threshold desired
beta = 0.8;
% Plotting of target tube
figure(1)
clf
hold on    
for time_indx = 0:time_horizon
    target_tube_at_time_indx = Polyhedron('H',...
        [target_tube(time_indx+1).A, ...
        zeros(size(target_tube(time_indx+1).A,1),1), ...
        target_tube(time_indx+1).b], 'He',[0 0 1 time_indx]);
    plot(target_tube_at_time_indx, 'alpha',0.25);
end
axis([-1 1 -1 1 0 time_horizon]);    
box on;
grid on;
xlabel('x');
ylabel('y');
zlabel('time');
title('Target tube');
axis equal;

%% Lagrangian underapproximation for stochastic reachability of a target tube
% |SReachSet| can compute Lagrangian underapproximation via multiple
% approaches. It can approximate the disturbance set, against which the robust
% computation is done, as a polytope or an ellipsoid (as specified by
% |bound_set_method|). Further, it can either perform an exact (and
% time-consuming) vertex-facet enumeration or an underapproximative (and
% complexity-preserving) computation using support vectors (as specified by
% |compute_style|). The vertex-facet enumeration may be performed by CDD or LRS,
% two popular enumeration techinques. We use CDD by default, but it can be
% switched to LRS using |vf_enum_method|.
n_dim = sys.state_dim + sys.input_dim; % Require unit vectors to sample X x U
lagrange_under_time = zeros(4,1);
%% Underapprox. type 1: Bound_set_method - Ellipsoid | Compute_style - VFmethod
% This is the best method in this case, since it leverages the fact that the
% distubance is Gaussian and for two-dimensional polytopes, vertex-facet
% enumeration is not significantly hard.
timerVal=tic;
luOpts1 = SReachSetOptions('term', 'lag-under', 'bound_set_method', ...
    'ellipsoid', 'verbose',1,'compute_style','vfmethod');
luOpts_time(1) = toc(timerVal);
timerVal=tic;
luSet(1) = SReachSet('term', 'lag-under', sys, beta, target_tube, luOpts1);
lagrange_under_time(1) = toc(timerVal);
%% Underapprox. type 2: Bound_set_method - Ellipsoid | Compute_style - Support
% Compared to type 1, this method will result in more conservative
% solution. The benefit of this approach becomes clear when we use it in
% high-dimensional systems where vertex-facet enumeration is hard. By the virtue
% of being recursion-free, this approach will scale better than type 1.
timerVal=tic;
luOpts2 = SReachSetOptions('term', 'lag-under', 'bound_set_method', ...
    'ellipsoid', 'system', sys, 'n_vertices', 2^n_dim*10+2*n_dim,...
    'verbose',1,'compute_style','support');
luOpts_time(2) = toc(timerVal);
timerVal=tic;
luSet(2) = SReachSet('term', 'lag-under', sys, beta, target_tube, luOpts2);
lagrange_under_time(2) = toc(timerVal);
%% Underapprox. type 3: Bound_set_method - Polytope | Compute_style - VFmethod
% Compared to type 1, this method will result in more conservative
% solution since the bounded disturbance set is larger in volume. This approach
% is best used when the disturbance is non-Gaussian. This approach is best used
% when the disturbance is non-Gaussian.
luOpts3 = SReachSetOptions('term', 'lag-under', 'bound_set_method', ...
    'polytope', 'verbose', 1, 'compute_style','vfmethod','template_polytope',...
    Polyhedron('lb',-ones(sys.dist.dim,1),'ub',ones(sys.dist.dim,1)));
luOpts_time(3) = toc(timerVal);
timerVal=tic;
luSet(3) = SReachSet('term', 'lag-under', sys, beta, target_tube, luOpts3);
lagrange_under_time(3) = toc(timerVal);
%% Underapprox. type 3: Bound_set_method - Polytope | Compute_style - Support
% Compared to type 1, this method will result in more conservative
% solution since the bounded disturbance set is larger in volume. Compared to
% type 3, this set will be more conservative since an underapproximative
% vertex-facet enumeration is performed.
luOpts4 = SReachSetOptions('term', 'lag-under', 'bound_set_method', ...
    'polytope', 'system', sys, 'n_vertices', 2^n_dim*10+2*n_dim,...
    'verbose', 1, 'compute_style','support', 'template_polytope',...
    Polyhedron('lb',-ones(sys.dist.dim,1),'ub',ones(sys.dist.dim,1)));
luOpts_time(4) = toc(timerVal);
timerVal=tic;
luSet(4) = SReachSet('term', 'lag-under', sys, beta, target_tube, luOpts4);
lagrange_under_time(4) = toc(timerVal);

%% Lagrangian overapproximation for stochastic reachability of a target tube
% |SReachSet| can compute Lagrangian overapproximation via multiple
% approaches. It can approximate the disturbance set which augments the
% input space as a polytope or an ellipsoid (as specified by
% |bound_set_method|). Further, it can perform a recursive computation
% using vertex-facet enumeration or a recursion-free computation using
% support functions (as specified by |compute_style|). The vertex-facet
% enumeration may be performed by CDD or LRS, two popular enumeration 
% techinques. We use CDD by default, but it can be switched to LRS using
% |vf_enum_method|.
n_dim_over = sys.state_dim;     % Require unit vectors to sample X
lagrange_over_time = zeros(4,1);
%% Overapprox. type 1: Bound_set_method - Ellipsoid | Compute_style - VFmethod
% This is the best method in this case, since it leverages the fact that the
% distubance is Gaussian and for two-dimensional polytopes, vertex-facet
% enumeration is not significantly hard.
timerVal=tic;
loOpts1 = SReachSetOptions('term', 'lag-over', 'bound_set_method', ...
    'ellipsoid', 'verbose', 1, 'compute_style','vfmethod');
loOpts_time(1) = toc(timerVal);
timerVal=tic;
loSet(1) = SReachSet('term', 'lag-over', sys, beta, target_tube, loOpts1);
lagrange_over_time(1) = toc(timerVal);
%% Overapprox. type 2: Bound_set_method - Ellipsoid | Compute_style - Support
% Compared to type 1, this method will provide an overapproximative solution but
% can potentially be faster since it is recursion-free.
timerVal=tic;
loOpts2 = SReachSetOptions('term', 'lag-over', 'bound_set_method', ...
    'ellipsoid', 'verbose', 1, 'compute_style','support', 'system', sys,...
    'n_vertices', 2^n_dim_over * 7+2*n_dim_over);
loOpts_time(2) = toc(timerVal);
timerVal=tic;
loSet(2) = SReachSet('term', 'lag-over', sys, beta, target_tube, loOpts2);
lagrange_over_time(2) = toc(timerVal);
%% Overapprox. type 3: Bound_set_method - Polytope | Compute_style - VFmethod
% Compared to type 1, this method will result in more conservative
% solution since the bounded disturbance set is larger in volume.
timerVal=tic;
loOpts3 = SReachSetOptions('term', 'lag-over', 'bound_set_method', ...
    'polytope', 'verbose', 1, 'template_polytope',...
    Polyhedron('lb',-ones(sys.dist.dim,1),'ub',ones(sys.dist.dim,1)),...
    'compute_style','vfmethod');
loOpts_time(3) = toc(timerVal);
timerVal=tic;
loSet(3) = SReachSet('term', 'lag-over', sys, beta, target_tube, loOpts3);
lagrange_over_time(3) = toc(timerVal);
%% Overapprox. type 4: Bound_set_method - Polytope | Compute_style - Support
% Compared to type 1, this method will result in more conservative
% solution since the bounded disturbance set is larger in volume. Compared to
% type 3, this method will provide an overapproximative solution but can
% potentially be faster since it is recursion-free.
loOpts4 = SReachSetOptions('term', 'lag-over', 'bound_set_method', ...
    'polytope', 'verbose', 1, 'template_polytope',...
    Polyhedron('lb',-ones(sys.dist.dim,1),'ub',ones(sys.dist.dim,1)),...
    'compute_style','support', 'system', sys,...
    'n_vertices', 2^n_dim_over* 7 +2*n_dim_over);
loOpts_time(4) = toc(timerVal);
timerVal=tic;
loSet(4) = SReachSet('term', 'lag-over', sys, beta, target_tube, loOpts4);
lagrange_over_time(4) = toc(timerVal);
                   

%% Dynamic programming solution
% We compare the results with <https://doi.org/10.1016/j.automatica.2010.08.006 
% dynamic programming> to see how the approximations appear and how they compare 
% in simulation times.
%%
dyn_prog_xinc = 0.025;
dyn_prog_uinc = 0.1;
tic;
[prob_x, cell_of_xvec] = SReachDynProg('term', sys, dyn_prog_xinc, ...
    dyn_prog_uinc, target_tube);
dynprog_time = toc();
%%
% Compute the beta-stochastic level set
%
dyn_soln_lvl_set=getDynProgLevelSets2D(cell_of_xvec, prob_x, beta, target_tube);
%% Simulation times: Lagrangian approximation beats dynamic programming 
% The simulation times for Lagrangian computation is much faster (in most
% cases) than dynamic programming. Further, dynamic programming solution
% provides no approximation guarantees while Lagrangian approach provides
% grid-free approximation guarantee.
%%
fprintf('Simulation times [seconds]:\n');
fprintf(' Lagrangian:\n');
% fprintf('   Overapproximation  : %.3f (online: %1.3f | offline: %1.3f)\n',...
%     lagrange_over_time + loOpts_time, lagrange_over_time, loOpts_time);
fprintf('   Overapproximation   online  |  offline | Total\n');
fprintf('(Ellipsoid,VFmethod)  %1.2e | %1.2e | %1.2f\n',...
    lagrange_over_time(1),loOpts_time(1),lagrange_over_time(1)+loOpts_time(1));
fprintf('(Ellipsoid, support)  %1.2e | %1.2e | %1.2f\n',...
    lagrange_over_time(2),loOpts_time(2),lagrange_over_time(2)+loOpts_time(2));
fprintf('(Polytope ,VFmethod)  %1.2e | %1.2e | %1.2f\n',...
    lagrange_over_time(3),loOpts_time(3),lagrange_over_time(3)+loOpts_time(3));
fprintf('(Polytope , support)  %1.2e | %1.2e | %1.2f\n',...
    lagrange_over_time(4),loOpts_time(4),lagrange_over_time(4)+loOpts_time(4));
fprintf('   Underapproximation  online  |  offline | Total\n');
fprintf('(Ellipsoid,VFmethod)  %1.2e | %1.2e | %1.2f\n',...
    lagrange_under_time(1),luOpts_time(1),lagrange_under_time(1)+luOpts_time(1));
fprintf('(Ellipsoid, support)  %1.2e | %1.2e | %1.2f\n',...
    lagrange_under_time(2),luOpts_time(2),lagrange_under_time(2)+luOpts_time(2));
fprintf('(Polytope ,VFmethod)  %1.2e | %1.2e | %1.2f\n',...
    lagrange_under_time(3),luOpts_time(3),lagrange_under_time(3)+luOpts_time(3));
fprintf('(Polytope , support)  %1.2e | %1.2e | %1.2f\n',...
    lagrange_under_time(4),luOpts_time(4),lagrange_under_time(4)+luOpts_time(4));
fprintf('   Dynamic programming: %.3f\n', dynprog_time);

%% Plotting all the sets together
% As expected, the over-approximation and the under-approximation obtained via 
% Lagrangian approach bounds the dynamic programming solution from "inside" and 
% "outside".
%%
figure(2);
clf
plot(safe_set, 'color', 'k');
hold on;
plot(loSet(2), 'color', 'y','alpha',1);
plot(loSet(1), 'color', 'r','alpha',0.5);
plot(dyn_soln_lvl_set,'color', 'b')
plot(luSet(1), 'color', 'm','alpha',1);
plot(luSet(2), 'color', 'g','alpha',0.75);
hold off;
xlabel('$x_1$', 'Interpreter', 'latex')
ylabel('$x_2$', 'Interpreter', 'latex')
leg = legend('Safe set',...
    'Overapproximation (Ell, support)',...
    'Overapproximation (Ell, VF)',...
     'Dyn. prog. soln.',...
    'Underapproximation (Ell, VF)',...
    'Underapproximation (Ell, support)');
set(leg,'Location','EastOutside');
box on;
axis equal;
axis tight;

figure(3);
clf
plot(safe_set, 'color', 'k');
hold on;
plot(loSet(4), 'color', 'y','alpha',1);
plot(loSet(3), 'color', 'r','alpha',0.5);
plot(dyn_soln_lvl_set,'color', 'b')
plot(luSet(3), 'color', 'm','alpha',1);
plot(luSet(4), 'color', 'g','alpha',0.75);
hold off;
xlabel('$x_1$', 'Interpreter', 'latex')
ylabel('$x_2$', 'Interpreter', 'latex')
leg = legend('Safe set',...
    'Overapproximation (Poly, support)',...
    'Overapproximation (Poly, VF)',...
     'Dyn. prog. soln.',...
    'Underapproximation (Poly, VF)',...
    'Underapproximation (Poly, support)');
set(leg,'Location','EastOutside');
box on;
axis equal;
axis tight;
