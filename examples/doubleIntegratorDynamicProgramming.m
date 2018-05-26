% SReachTools/doubleIntegratorDynamicProgramming
% ============================================================================
% 
% Dyhamic programming for reachability of a target tube.
%
% This example will compute the dynamic programming solution for staying in a 
% safe set `safe_set = Polyhedron('lb', [-1, -1], 'ub', [1, 1])` (the viability
% problem). The dynamic are discretized double integrator dynamics:
%
%   x_{k+1} = [ 1, T;  x_{k} + [ T^2/2;   u_{k} + w_{k}
%               0, 1];           T     ]; 
% 
% w_{k} is a Guassian diturbance, w_{k} ~ N([0, 0]', I_{2})
%
% ============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/abyvinod/SReachTools/blob/master/LICENSE
% 
% 

% example parameters
T = 0.25;

% define the system
sys = LtiSystem('StateMatrix', [1, T; 0, 1], ...
    'InputMatrix', [T^2/2; T], ...
    'InputSpace', Polyhedron('lb', -0.1, 'ub', 0.1), ...
    'DisturbanceMatrix', eye(2), ...
    'Disturbance', StochasticDisturbance('Gaussian', zeros(2,1), 0.01*eye(2)));

% solving the "viability" problem means we want the system to stay within a set
% of safe states
% safe set definition
safe_set = Polyhedron('lb', [-1, -1], 'ub', [1, 1]);
target_set = Polyhedron('lb', [-0.5, -0.5], 'ub', [0.5, 0.5]);

% in target tube for the viability problem is equivalent to a tube of repeating
% safe sets
target_tube = {safe_set, ...
    safe_set, ...
    safe_set, ...
    safe_set, ...
    safe_set, ...
    safe_set, ...
    safe_set, ...
    safe_set, ...
    safe_set, ...
    safe_set, ...
    target_set};

N = length(target_tube);

% need to create a state space grid and input space grid
ss_grid = SpaceGrid([-1, -1], [1, 1], 40);
in_grid = InputGrid(-1, 1, 20);

grid_probability = getDynProgSolForTargetTube(sys, ...
    ss_grid, in_grid, target_tube, 'DisplayStatus', 3);

%% plotting
figure(1);
ss_grid.plotGridProbability(grid_probability);
view(0, 90)


%% Moving target problem
target_tube = {Polyhedron('lb', [-1, -1], 'ub', [1, 1]), ...
    Polyhedron('lb', [-0.5, -1], 'ub', [1, 0.5]), ...
    Polyhedron('lb', [-1, -1], 'ub', [0.5, 0.5]), ...
    Polyhedron('lb', [-1, -0.5], 'ub', [0.5, 1]), ...
    Polyhedron('lb', [-0.5, -0.5], 'ub', [1, 1])};

ss_grid = SpaceGrid([-1, -1], [1, 1], 40);
in_grid = InputGrid(-0.1, 0.1, 3);

grid_probability = getDynProgSolForTargetTube(sys, ...
    ss_grid, in_grid, target_tube);

figure(2);
ss_grid.plotGridProbability(grid_probability);
view(0, 90)