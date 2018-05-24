function sys = getChainOfIntegLtiSystem(dim, T, inputSpace, dist_mat, distSpace)
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
%       'InputSpace', Polyhedron('lb', -1, 'ub', 1));
%
% ============================================================================
% 
% sys = getChainOfIntegLtiSystem(dim, T, input_space, disturb)
% 
% Inputs:
% -------
%   dim         - Dimensions
%   T           - Discretization time step
%   input_space - Input space (Polyhedron)
%   disturb     - Disturbance object (Polyhedron / StochasticDisturbance)
% 
% Outputs:
% --------
%   sys - LtiSystem object
% 
% =============================================================================
%
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 
%

    inpar = inputParser();
    inpar.addRequired('dim', @(x) validateattributes(x, ...
        {'numeric'}, {'integer', '>', 0}));
    inpar.addRequired('T', @(x) validateattributes(x, ...
        {'numeric'}, {'>', 0}));
    inpar.addRequired('inputSpace', @(x) validateattributes(x,...
        {'Polyhedron'}, {'nonempty'}));
    inpar.addRequired('dist_mat', @(x) validateattributes(x, {'numeric'},...
        {'nonempty'}));
    inpar.addRequired('distSpace',  @(x) validateattributes(x,...
        {'Polyhedron', 'StochasticDisturbance', 'RandomVector'},{'nonempty'}));

    try
        inpar.parse(dim, T, inputSpace, dist_mat, distSpace)
    catch err
        exc = SrtInvalidArgsError.withFunctionName();
        exc = exc.addCause(err);
        throwAsCaller(exc);
    end

    % if input space is defined it can only be of dimension 1
    if inpar.Results.inputSpace.Dim > 1
        exc = SrtInvalidArgsError();
        exc = exc.addCause(SrtInvalidArgsError(['Chain of integrator ', ...
            'input space can only have dimension of 1']));
        throwAsCaller(exc);
    end

    if xor(any(strcmp('Disturbance', inpar.UsingDefaults)), ...
           any(strcmp('DisturbanceMatrix', inpar.UsingDefaults)))
        exc = SrtInvalidArgsError();
        exc = exc.addCause(SrtInvalidArgsError(['Disturbance ', ...
            'matrix and disturbance must be simultaneously defined']));
        throwAsCaller(exc);
    end

    % anonymous function for getting the necessary matrix internals
    facT = @(t, n) t^n / factorial(n);

    % initialization
    state_mat = eye(dim);
    input_mat = zeros(dim, 1);
    
    % Populate the upper triangle of state_mat and the entries of input_mat
    for i = 1:dim
        input_mat(i) = facT(T, dim-i+1);
        for j = i+1:dim
            state_mat(i, j) = facT(T, j-i);
        end
    end

    sys = LtiSystem('StateMatrix', state_mat, ...
        'InputMatrix', input_mat, ...
        'InputSpace', inpar.Results.inputSpace, ...
        'DisturbanceMatrix', inpar.Results.dist_mat, ...
        'Disturbance', inpar.Results.distSpace);
end
