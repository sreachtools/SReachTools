function sys = getChainOfIntegLtiSystem(dim, T, input_space, disturb)
% SReachTools/systems/getChainofIntegLtiSystem  Get chain of integrators Lti 
% System
% ============================================================================
% 
% Get an n-d discrete chain of integrators in the form of an LtiSystem object
% given a discretization time-step T.
%
% Chain of integrators model:
%
%   x_{k+1} = A * x_{k} + B * u_{k} + w_{k}
% 
%   A = [1, T, T^2/2, ... T^n/n!;
%        0, 1, T,     ... T^(n-1)/(n-1)!;
%        ...
%        0, 0, ...    ... T;
%        0, 0, ...    ... 1];
%  
%   B = [T^n/n!;
%        T^(n-1)/(n-1)!
%        ...
%        T];
% 
% Usage:
% ------
%   % 3-d chain of integrators with U = [-1,1] and no (empty) disturbance
%   sys = getChainOfIntegLtiSystem(3, 0.2, ...
%       Polyhedron('lb', -1, 'ub', 1), ...
%       Polyhedron());
%
% ============================================================================
% 
% sys = getChainOfIntegLtiSystem(dim, T, input_space, disturb)
% 
% Inputs:
% -------
%   dim         - Dimensions
%   T           - Discretization time step
%   input_space - Input space
%   disturb     - Disturbance
% 
% Outputs:
% --------
%   sys - LtiSystem object
% 
% =============================================================================
%
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/abyvinod/SReachTools/blob/master/LICENSE
% 

% anonymous function for getting the necessary matrix internals
facT = @(t, n) t^n / factorial(n);

% initialization
state_mat = eye(dim);
input_mat = zeros(dim, 1);

for i = 1:dim
    input_mat(i) = facT(T, dim-i+1);
    for j = i+1:dim
        state_mat(i, j) = facT(T, j-i);
    end
end

sys = LtiSystem('StateMatrix', state_mat, ...
    'InputMatrix', input_mat, ...
    'InputSpace', input_space, ...
    'DisturbanceMatrix', eye(2), ...
    'Disturbance', disturb);

end