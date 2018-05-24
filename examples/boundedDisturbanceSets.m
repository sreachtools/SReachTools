% SReachTools/boundedDisturbanceSets: Bounded disturbance set generation example
% ============================================================================
%
% Example script for generating bounded disturbances with different methods.
% Specific mehods available:
%   random    - Ellipsoid beneration with randomly chosen directions
%   box       - N-d box
% 
% ============================================================================
% 
% This function is part of the Stochastic Optimal Control Toolbox.
% License for the use of this function is given in
%      https://github.com/abyvinod/SReachTools/blob/master/LICENSE
% 
% 

d = StochasticDisturbance('Gaussian', zeros(2,1), 0.005*eye(2));

p1 = getBoundedSetForDisturbance(d, 3, 0.8, 'random', 50);
p2 = getBoundedSetForDisturbance(d, 3, 0.8, 'box', 1e-4);

figure(1)
plot(p1, 'Color', 'y')
hold on;
plot(p2)
hold off;
