% SReachTools/doubleIntegratorMovingTarget
% ============================================================================
% 
% This script is an example application of using SReachTools to perform a dynamic 
% programming solution to a reachability of a target tube problem and using
% the Lagrangian under and overapproximation methods. The example uses double 
% integrator dynamics:
%   
%   x_{k+1} = [ 1, T;  x_{k} + [ T^2/2;   u_{k} + w_{k}
%               0, 1];   
%
% With a Gaussian disturbance w_{k} ~ N([0, 0]', 0.002 * I).
% 
% The target tube rotates around the box |x| < 1, |y| < 1. The order of rotation
% is NE, NW, SW, SE, entire region.
%
% ============================================================================
% 
% This function is part of the Stochastic Optimal Control Toolbox.
% License for the use of this function is given in
%      https://github.com/abyvinod/SReachTools/blob/master/LICENSE
% 
%

% example parameters
T = 0.25;

% define the system
sys = getChainOfIntegLtiSystem(2, ...
    T, ...
    Polyhedron('lb', -0.1, 'ub', 0.1), ...
    StochasticDisturbance('Gaussian', zeros(2,1), 0.002*eye(2)));

%% Moving target problem
target_tube = {Polyhedron('lb', [-1, -1], 'ub', [1, 1]), ...
    Polyhedron('lb', [-0.5, -1], 'ub', [1, 0.5]),...
    Polyhedron('lb', [-1, -1], 'ub', [0.5, 0.5]), ...
    Polyhedron('lb', [-1, -0.5], 'ub', [0.5, 1]), ...
    Polyhedron('lb', [-0.5, -0.5], 'ub', [1, 1])};

ss_grid = SpaceGrid([-1, -1], [1, 1], 100);
in_grid = InputGrid(-0.1, 0.1, 3);

n_targets = length(target_tube);

% dynamic programming
figure()
[X, Y] = ss_grid.getMeshGrids();
for i = 1:n_targets
    grid_probability = getDynProgSolForTargetTube(sys, ...
        ss_grid, in_grid, target_tube(n_targets-i+1:n_targets));
    
    under_approx_set = getApproxStochasticLevelSetViaLagrangian(sys, ...
        0.8, ...
        target_tube(n_targets-i+1:n_targets), ...
        'underapproximation', ...
        'random', 50);
    
    over_approx_set = getApproxStochasticLevelSetViaLagrangian(sys, ...
        0.8, ...
        target_tube(n_targets-i+1:n_targets), ...
        'overapproximation', ...
        'random', 50);
    
    subplot(1, n_targets, i)
    plot(target_tube{n_targets-i+1}, 'Color', 'y');
    hold on;
    plot(over_approx_set, 'Color', 'r');
    contourf(X, Y, reshape(grid_probability, ss_grid.n_points), ...
        'LevelList', 0.8)
    caxis([0.8, 1])
    plot(under_approx_set, 'Color', 'g')
    hold off;
end
