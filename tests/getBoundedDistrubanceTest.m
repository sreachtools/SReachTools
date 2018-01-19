% Description : Unit test script for the function
%               getBoundedSetForDisturbance
% 1/30/2018
%

%% Bounded Disturbance Ellipse for Gaussian Disturbance with Random Directions
% double integrator system
time_horizon = 5;
T = 0.25;
safe_max = 1;
target_max = 0.5;
umax = 0.75;
dmax = 0.097;

sys = LtiSystem(...
    'StateMatrix', [1, T; 0, 1], ...
    'InputMatrix', [T^2; T], ...
    'InputSpace', Polyhedron('lb', -umax, 'ub', umax), ...
    'DisturbanceMatrix', eye(2), ...
    'Disturbance', StochasticDisturbance('Gaussian', zeros(2,1), 0.01* eye(2)));

beta = 0.7;
had_error = 0;
try
    % get bounded set
    bounded_set = getBoundedSetForDisturbance(...
        sys.disturbance, ...
        time_horizon, ...
        beta, ...
        'random', ...
        100); 
    
    % since the choice of directions is random we can't compare with an existing
    % polytope so best I can do is try to validate that the result is a
    % polyhedron
    validateattributes(bounded_set, {'Polyhedron'}, {'nonempty'})
catch ME
    had_error = 1;
    throw(ME)
end

assert(~had_error, ['Error occurred incomputing the bounding set with ', ...
    'random method']);

%% Load Bounded Disturbance Ellipse for CWH Problem
% load cwh distrubance ellipsoid

successful = true;
try
    % get bounded set
    bounded_set = getBoundedSetForDisturbance(...
        [], ...
        [], ...
        [], ...
        'load', ...
        'data/cwhUnderapproxBoundeSet.mat'); 
    
    % since the choice of directions is random we can't compare with an existing
    % polytope so best I can do is try to validate that the result is a
    % polyhedron
    validateattributes(bounded_set, {'Polyhedron'}, {'nonempty'})
catch ME
    successful = false;
    throw(ME)
end

assert(successful, 'Error occurred in loadingbounding set');