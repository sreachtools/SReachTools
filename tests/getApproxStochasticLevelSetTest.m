% Description : Unit test script for the function
%               getApproxStochasticLevelSEtViaLagrangian
% 1/30/2018
%

%% Underapproximation of double integrator via lagrangian
% double integrator dynamics
T = 0.25;

% Safe set K |x1| < 1, |x2| < 1
safe_set = Polyhedron('lb', [-1; -1], 'ub', [1; 1]);
time_horizon = 6;

% Input Space
U = Polyhedron('lb', -0.1, 'ub', 0.1);

sys = LtiSystem('StateMatrix', [1, T; 0, 1], ...
    'InputMatrix', [T^2/2; T], ...
    'InputSpace', U, ...
    'DisturbanceMatrix', eye(2), ...
    'Disturbance', RandomVector('Gaussian', zeros(2,1), 5e-3*eye(2)));

% target_tube = {K, K, K, K, K, K};
target_tube = TargetTube('viability', safe_set, time_horizon);

successful = true;
try
    level_set = getApproxStochasticLevelSetViaLagrangian(sys, 0.8, ...
        target_tube, 'underapproximation', 'random', 100);
catch ME
    successful = false;
    throw(ME)
end

assert(successful, ['Error in computing underapproximation of stochastic ', ...
    'level set for double integrator'])