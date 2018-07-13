function sys = getCwhLtiSystem(dim, varargin)
% SReachTools/systemDefinitions/getCwhLtiSystem: Create a LtiSystem object for
% the spacecraft dynamics using Clohessy-Wiltshire-Hill (CWH) dynamics
% =============================================================================
%
% Constructs a LtiSystem object for the discrete-time linear time-invariant
% dynamics of an approaching spacecraft (the deputy) relative to a target (the
% chief).
% 
% In addition, you can provide:
% - an input space for the control inputs are the components of the external
%   force vector, i.e. the thruster control input.  
% - an additive Gaussian noise process noise vector that represents the
%   uncertainty in the model due to external forces on the spacecraft not
%   captured in the linearized model.  
% - Other parameters relevant to the Clohessy-Wiltshire-Hill (CWH) dynamics in
%   the form of a struct (see Notes) 
%
% The continuous-time is given in (1),(2) of
%      K. Lesser, M. Oishi, and R. S. Erwin, "Stochastic Reachability for
%      Control of Spacecraft Relative Motion", in Proceedings of IEEE
%      Conference on Decision and Control, 2013. 
% Alternatively, see 
%   - Curtis, Howard D. Orbital mechanics for engineering students.
%     Butterworth-Heinemann, 2013, Section 7.4
%
% The state of the system is [position in x, position in y, position in z,
% velocity in x, velocity in y, and velocity in z].
%
% Usage:
% ------
%
% % Create a LtiSystem for the CWH dynamics using the parameters given in Lesser
% % et. al, CDC 2013 paper.
%
% sys = getCwhLtiSystem(4, ...
%                       Polyhedron('lb', -0.01*ones(2,1), ...
%                                  'ub',  0.01*ones(2,1)), ...
%                       StochasticDisturbance(...
%                             'Gaussian', ...
%                             zeros(4,1), ...
%                             diag([1e-4, 1e-4, 5e-8, 5e-8])));
%
% % Create a LtiSystem for the uncontrolled 6D CWH dynamics
% sys = getCwhLtiSystem(6);
%
% % Create a LtiSystem for the uncontrolled 4D CWH dynamics
% sys = getCwhLtiSystem(4, Polyhedron(), StochasticDisturbance(...
%                                             'Gaussian', ...
%                                             zeros(4,1), ...
%                                             diag([1e-4, 1e-4, 5e-8, 5e-8])));
%
% % Create a LtiSystem for the controlled 6D CWH dynamics
% sys = getCwhLtiSystem(6, ...
%                       Polyhedron('lb', -0.01*ones(3,1), ...
%                                  'ub',  0.01*ones(3,1)), ...
%                       StochasticDisturbance(...
%                             'Gaussian', ...
%                             zeros(6,1), ...
%                             diag([1e-4, 1e-4, 1e-4, 5e-8, 5e-8, 5e-8])));
%
% =============================================================================
%
% Inputs:
% -------
%   dim         - Dimension of the CWH dynamics of interest (Needs to be 4 or 6)
%   input_space - (Optional) Input space for the spacecraft (Polytope)
%                 [Provide an empty polyhedron to create an uncontrolled but
%                  perturbed system]
%   disturbance - (Optional) Stochastic disturbance object describing the
%                 disturbance affecting the dynamics
%   user_params - (Optional) User parameter struct that gives as a name-value
%                 pair different parameters affecting the dynamics.  Possible
%                 values that may be adjusted are --- sampling_period,
%                 orbital_radius, grav_constant, celes_mass, chief_mass,
%                 orbit_ang_vel, disc_orbit_dist 
%                 [If empty, default values are set.]
%
% Outputs:
% --------
%   sys - LtiSystem object describing the CWH dynamics
%
% Notes:
% ------
% * This code and the parameters were obtained from Lesser's repeatability code
%   for the 2013 CDC paper.
% * The default parameters for the CWH system dynamics are:
%       sampling period              = 20 s
%       orbital radius               = 850 + 6378.1 m
%       gravitational constant       = 6.673e-11
%       celestial body mass          = 5.9472e24 kg
%       gravitational body           = grav_constant * celes_mass / 1e6
%       orbital angular velocity     = sqrt(grav_body / orbital_radius^3)
%       chief mass                   = 300 kg
%       discretized orbital distance = orbit_ang_vel * sampling_period rad
%
% =============================================================================
%
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 

    %% CWH dynamics parameters --- Default or user-provided
    if length(varargin) == 3
        params = checkSystemParameters(varargin{3});
    elseif length(varargin) < 3
        params = getDefaultCwhParameters();
    else
        error('SReachTools:invalidArgs', 'Too many inputs');
    end
    
    % System matrix construction
    [state_matrix_6D, input_matrix_6D] = get6dCwhStateAndInputMatrices(params);
    if dim == 4
        state_matrix = state_matrix_6D([1:2,4:5],[1:2,4:5]);
        input_matrix = input_matrix_6D([1:2,4:5],1:2);
    elseif dim == 6
        state_matrix = state_matrix_6D;
        input_matrix = input_matrix_6D;
    else
        error('SReachTools:invalidArgs', 'dim must be 4 or 6.');
    end
    
    
    %% Input handling
    if isempty(varargin)
        % Construct a LtiSystem object without disturbance or input space
        sys = LtiSystem('StateMatrix', state_matrix);
    elseif length(varargin) == 1
        % Only input space --- No disturbance
        assert( isa(varargin{1}, 'Polyhedron'), ...
               'SReachTools:invalidArgs', ...
               'Must provide polyhedral input space');           
        % Construct a LtiSystem object without disturbance
        sys = LtiSystem('StateMatrix', state_matrix, ...
                        'InputMatrix', input_matrix, ...
                        'InputSpace', varargin{1});
    else        % Since length(varargin) < 3, length(varargin) == 2
        % Input space and disturbance
        assert( isa(varargin{1}, 'Polyhedron'), ...
               'SReachTools:invalidArgs', ...
               'Must provide polyhedral input space');           
        assert( isa(varargin{2}, 'RandomVector'), ...
                'SReachTools:invalidArgs', ...
                'Must provide a random vector object');
        % Disturbance matrix construction is nxn identity
        disturbance_matrix = eye(size(state_matrix,2));
        if isEmptySet(varargin{1})
            % Construct a LtiSystem object without input space
            sys = LtiSystem('StateMatrix', state_matrix, ...
                            'DisturbanceMatrix', disturbance_matrix, ...
                            'Disturbance', varargin{2});            
        else
            % Construct a LtiSystem object with everything in
            sys = LtiSystem('StateMatrix', state_matrix, ...
                            'InputMatrix', input_matrix, ...
                            'InputSpace', varargin{1}, ...
                            'DisturbanceMatrix', disturbance_matrix, ...
                            'Disturbance', varargin{2});            
        end
    end    
end

function [state_matrix, input_matrix] = get6dCwhStateAndInputMatrices(params)
% SReachTools/getCwhLtiSystem/get6dCwhStateAndInputMatrices: Get 6-d CWH 
% matrices for discrete-time LTI CWH dynamics
% =============================================================================
%
% Get the discretized state and input matrices for the 6-d CWH system
%   
% Usage: Nested function
%
% =============================================================================
%
% Inputs:
% -------
%   params - Parameter struct
%
% Outputs:
% --------
%   state_matrix - Discrete-time LTI CWH state matrix
%   input_matrix - Discrete-time LTI CWH input matrix
%
% =============================================================================
%
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 

    % redefine important variable for simplicity
    sampling_period = params.sampling_period;
    orbit_ang_vel = params.orbit_ang_vel;
    chief_mass = params.chief_mass;

    % Continuous-time LTI CWH unforced dynamics e^{A_{cts}t}
    e_power_A_cts_t = @(t) [ 
        4 - 3 * cos(orbit_ang_vel * t), ...
            0, ...
            0, ...
            (1/orbit_ang_vel) * sin(orbit_ang_vel * t), ...
            (2/orbit_ang_vel) * (1 - cos(orbit_ang_vel * t)), ...
            0; ...
        6 * (sin(orbit_ang_vel * t) - orbit_ang_vel * t), ...
            1, ...
            0, ...
            -(2/orbit_ang_vel) * (1 - cos(orbit_ang_vel * t)), ...
            (1/orbit_ang_vel) * (4*sin(orbit_ang_vel * t) - 3*orbit_ang_vel * t), ...
            0; ...
        0, ...
            0, ...
            cos(orbit_ang_vel * t), ...
            0, ...
            0, ...
            (1/orbit_ang_vel) * sin(orbit_ang_vel * t); ...
        3 * orbit_ang_vel * sin(orbit_ang_vel * t), ...
            0, ...
            0, ...
            cos(orbit_ang_vel * t), ...
            2 * sin(orbit_ang_vel * t), ...
            0; ...
        -6 * orbit_ang_vel * (1 - cos(orbit_ang_vel * t)), ...
            0, ...
            0, ...
            -2 * sin(orbit_ang_vel * t), ...
            4 * cos(orbit_ang_vel * t) - 3, ...
            0; ...
        0, ...
            0, ...
            -orbit_ang_vel * sin(orbit_ang_vel * t), ...
            0, ...
            0, ...
            cos(orbit_ang_vel * t);
    ];
    % Discrete-time state matrix is Phi(T_s) for sampling time T_s since the
    % system is time-invariant
    state_matrix = e_power_A_cts_t(sampling_period);

    % Continuous-time input matrix B_{cts}
    B_cts = 1/chief_mass*[zeros(3);eye(3)];
    % Discrete-time input matrix is (\int_0^T e^{A_{cts}\tau} d\tau) B_cts
    input_matrix = integral(e_power_A_cts_t, ...
                            0, ...
                            sampling_period, ...
                            'ArrayValued', true) * B_cts;
end


function params = getDefaultCwhParameters()
% SReachTools/getCwhLtiSystem/getDefaultCwhParameters: Get default parameter
% struct
% =============================================================================
%
% Get the default parameters for the CWH system dynamics. Defaults:
%   sampling period              = 20 s
%   orbital radius               = 850 + 6378.1 m
%   gravitational constant       = 6.673e-11
%   celestial body mass          = 5.9472e24 kg
%   gravitational body           = grav_constant * celes_mass / 1e6
%   orbital angular velocity     = sqrt(grav_body / orbital_radius^3)
%   chief mass                   = 300 kg
%   discretized orbital distance = orbit_ang_vel * sampling_period rad
%   
% Usage: Nested function
%
% =============================================================================
%
% Inputs: None
%
% Outputs:
% --------
%   params - Parameter struct
%
% =============================================================================
%
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 

    params = struct();

    % sampling period in sec
    params.sampling_period = 20;

    % Orbital radius in m
    params.orbital_radius = 850 + 6378.1;

    % Universal gravitational constant                                           
    params.grav_constant = 6.673e-11;

    % Mass of the celestial body (default is Earth)
    params.celes_mass = 5.9472e24;

    % Gravitation constant for the pull of the celestial body (default Earth)
    params.grav_body = params.grav_constant * params.celes_mass / 1000^3;

    % Angular velocity in the orbit
    params.orbit_ang_vel = sqrt(params.grav_body / params.orbital_radius^3);

    % Mass of the chief kg
    params.chief_mass = 300; 

    % Discretized orbital distance
    params.disc_orbit_dist = params.orbit_ang_vel * params.sampling_period;                


end

function params = checkSystemParameters(user_params)
% SReachTools/getCwhLtiSystem/checkSystemParameters: Check the given system
% parameters
% =============================================================================
%
% Check the provided system parameters and update any defaults
%   
% Usage: Nested function
%
% =============================================================================
%
% Inputs:
% -------
%   user_params - User parameter struct
%                 Possible values that may be adjusted are --- sampling_period,
%                 orbital_radius, grav_constant, celes_mass, chief_mass,
%                 orbit_ang_vel, disc_orbit_dist
%                 If empty, default values are set.
%
% Outputs:
% --------
%   params - Parameter struct
%
% =============================================================================
%
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 

    params = getDefaultCwhParameters();
    all_fields = fields(params);

    user_fields = fields(params);
    for i = 1:length(all_fields)
        switch(lower(all_fields{i}))
            case 'sampling_period'
                if ismember('sampling_period', user_fields)
                    params.sampling_period = user_params.sampling_period;
                end
            case 'orbital_radius'
                if ismember('orbital_radius', user_fields)
                    params.orbital_radius = user_params.orbital_radius;
                end
            case 'grav_constant'
                if ismember('grav_constant', user_fields)
                    warning(['Gravitational constant is a universal ', ...
                        'constant, simulation for different universes not ', ...
                        'available. Using this universe''s constant of ', ...
                        '%f.'], params.grav_constant)
                end
            case 'celes_mass'
                if ismember('celes_mass', user_fields)
                    params.celes_mass = user_params.celes_mass;
                end
            case 'grav_body'
                if ismember('grav_body', user_fields)
                    warning(['The gravitational force of a celestial body ', ...
                        'is determined via an algebraic equation. Cannot ', ...
                        'set this value.'])
                end

                params.grav_body = params.grav_constant * params.celes_mass/...
                    1000^3;
            case 'orbit_ang_vel'
                if ismember('orbit_ang_vel', user_fields)
                    warning(['The orbital angular velocity of a satellite ', ...
                        'is determined via an algebraic equation. Cannot ', ...
                        'set this value.'])
                end

                params.orbit_ang_vel = sqrt(poarams.grav_body / ...
                    params.orbital_radius^3);
            case 'chief_mass'
                if ismember('chief_mass', user_fields)
                    params.chief_mass = user_params.chief_mass;
                end
            case 'disc_orbit_dist'
                if ismember('sampling_period', user_fields)
                    warning(['The discretized orbital distance of a ', ...
                        'satellite is determined via an algebraic ', ...
                        'equation. Cannot set this value.'])
                end

                params.disc_orbit_dist = params.orbit_ang_vel * ...
                    params.sampling_period;

            otherwise
                error('SReachTools:invalidArgs', 'Unhandled option');
        end
    end
end
