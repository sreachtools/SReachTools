function sys = getCwhLtiSystem(dim, ...
                               input_space, ...
                               disturbance, ...
                               user_params)
% SReachTools/systemDefinitions/getCwhLtiSystem: Create a LtiSystem object for
% the spacecraft dynamics using Clohessy-Wiltshire-Hill (CWH) dynamics
% =============================================================================
%
% Constructs a LtiSystem object for the in-plane (four-dimensional state),
% discretized in time, linearized, and time-invariant dynamics of an
% approaching spacecraft (the deputy) relative to a target (the chief). The
% control inputs are the components of the external force vector, i.e. the
% thruster control input. The model has an additive Gaussian noise process
% noise vector that represents the uncertainty in the model due to external
% forces on the spacecraft not captured in the linearized model.  The
% out-of-plane motion (not modeled in this LtiSystem object) is a
% two-dimensional problem (decoupled from this system) within the Hill
% reference frame.
%
% The continuous-time (which is later discretized in time) in-plane relative
% dynamics, known as Clohessy-Wiltshire-Hill equations, is given in (1),(2) of
%      K. Lesser, M. Oishi, and R. S. Erwin, "Stochastic Reachability for
%      Control of Spacecraft Relative Motion", in Proceedings of IEEE
%      Conference on Decision and Control, 2013. 
% Alternatively, see W. Wiesel, Spaceflight Dynamics. New York: McGraw-Hill,
% 1989.
%
% The state of the system is [position in x, position in y, velocity in x and
% velocity in y].
%
% Usage:
% ------
%
% % Create a LtiSystem for the CWH dynamics using the parameters given in Lesser
% % et. al, CDC 2013 paper.
%
% sys = getCwhLtiSystem(4,...
%                       Polyhedron('lb', -0.01*ones(2,1),...
%                                  'ub',  0.01*ones(2,1)),...
%                       StochasticDisturbance(...
%                             'Gaussian',...
%                             zeros(4,1),...
%                             diag([1e-4, 1e-4, 5e-8, 5e-8])));
%
% =============================================================================
%
% Inputs:
% -------
%   dim         - Dimension of the CWH dynamics of interest (4 or 6)
%   input_space - Input space for the spacecraft (Polytope)
%   disturbance - Stochastic disturbance object describing the disturbance
%                     affecting the dynamics
%   user_params - User parameter struct that gives as a name-value pair
%                 different parameters affecting the dynamics.
%                 Possible values that may be adjusted are --- sampling_period,
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
%   sampling period              = 20
%   orbital radius               = 850 + 6378.1
%   gravitational constant       = 6.673e-11
%   celestial body mass          = 5.9472e24
%   gravitational body           = grav_constant * celes_mass / 1000^3
%   orbital angular velolicty    = sqrt(grav_body / orbital_radius^3)
%   chief mass                   = 300 
%   discretized orbital distance = orbit_ang_vel * sampling_period
%
% =============================================================================
%
%   This function is part of the Stochastic Optimal Control Toolbox.
%   License for the use of this function is given in
%        https://github.com/abyvinod/SReachTools/blob/master/LICENSE
% 

    %% Input handling
    % Input must be a 2-dimensional polytope
    assert( isa(input_space, 'Polyhedron'), ...
           'SReachTools:invalidArgs',...
           'Must provide polyhedral input space');
    assert( isa(disturbance, 'StochasticDisturbance'), ...
            'SReachTools:invalidArgs', ...
            'Must provide a stochastic disturbance');

    %% CWH dynamics parameters
    if nargin < 3
        error('SReachTools:internal', 'Too few input arguments.');
    elseif nargin < 4
        params = getDefaultCwhParameters();
    else
        params = checkSystemParameters(user_params);
    end
    
    % System matrix construction
    if dim == 4
        [state_matrix, input_matrix] = get4dCwhStateAndInputMatrices(params);
    else
        [state_matrix, input_matrix] = get6dCwhStateAndInputMatrices(params);
    end

    % Disturbance matrix construction
    disturbance_matrix = eye(size(state_matrix,2));

    % Construct a LtiSystem object
    sys = LtiSystem('StateMatrix', state_matrix,...
                    'InputMatrix', input_matrix,...
                    'InputSpace', input_space,...
                    'DisturbanceMatrix', disturbance_matrix,...
                    'Disturbance', disturbance);
end

function params = getDefaultCwhParameters()
% SReachTools/getCwhLtiSystem/getDefaultCwhParameters: Get default parameter
% struct
% =============================================================================
%
% Get the default parameters for the CWH system dynamics. Defaults:
%   sampling period              = 20
%   orbital radius               = 850 + 6378.1
%   gravitational constant       = 6.673e-11
%   celestial body mass          = 5.9472e24
%   gravitational body           = grav_constant * celes_mass / 1000^3
%   orbital angular velolicty    = sqrt(grav_body / orbital_radius^3)
%   chief mass                   = 300 
%   discretized orbital distance = orbit_ang_vel * sampling_period
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
%   This function is part of the Stochastic Optimal Control Toolbox.
%   License for the use of this function is given in
%        https://github.com/abyvinod/SReachTools/blob/master/LICENSE
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
%   This function is part of the Stochastic Optimal Control Toolbox.
%   License for the use of this function is given in
%        https://github.com/abyvinod/SReachTools/blob/master/LICENSE
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

                params.disc_orbit_dist = params.orbit_ang_vel * params.sampling_period;

            otherwise
                assert(false, 'Unhandled option')
        end
    end
end

function [state_matrix, input_matrix] = get4dCwhStateAndInputMatrices(params)
% SReachTools/getCwhLtiSystem/get4dCwhStateAndInputMatrices: Get 4-d CWH 
% matrices
% =============================================================================
%
% Get the discretized state and input matrices for the 4-d CWH system
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
%   state_matrix - CWH state matrix
%   input_matrix - CWH input matrix
%
% =============================================================================
%
%   This function is part of the Stochastic Optimal Control Toolbox.
%   License for the use of this function is given in
%        https://github.com/abyvinod/SReachTools/blob/master/LICENSE
% 

    % redefine important variable for simplicity
    sampling_period = params.sampling_period;
    disc_orbit_dist = params.disc_orbit_dist;
    orbit_ang_vel = params.orbit_ang_vel;
    chief_mass = params.chief_mass;

    % System matrix construction
    state_matrix = [
        4 - 3*cos(disc_orbit_dist), ...
            0, ...
            sin(disc_orbit_dist)/orbit_ang_vel, ...
            (2 / orbit_ang_vel)*(1 - cos(disc_orbit_dist));
        6 * (sin(disc_orbit_dist) - disc_orbit_dist), ...
            1, ...
            - (2 / orbit_ang_vel) * (1-cos(disc_orbit_dist)),...
            (4 * sin(disc_orbit_dist) - 3 * disc_orbit_dist) / ...
                orbit_ang_vel;
        3 * orbit_ang_vel * sin(disc_orbit_dist), ...
            0, ...
            cos(disc_orbit_dist), ...
            2*sin(disc_orbit_dist);
        -6 * orbit_ang_vel * (1-cos(disc_orbit_dist)), ...
            0, ...
            -2 * sin(disc_orbit_dist), ...
            4 * cos(disc_orbit_dist) - 3
    ];


    % Input matrix construction
    int_eat = @(t) [
        4 * t - 3 * sin(orbit_ang_vel*t) / orbit_ang_vel, ...
            0,...
            -cos(orbit_ang_vel * t) / orbit_ang_vel^2, ...
            (2 / orbit_ang_vel) * ...
                (t - sin(orbit_ang_vel * t) / orbit_ang_vel);
        6 * (sin(orbit_ang_vel * t) / orbit_ang_vel - ...
            orbit_ang_vel*t^2/2), ...
            t,... 
            -(2 / orbit_ang_vel) * ...
                (t - sin(orbit_ang_vel * t) / orbit_ang_vel),... 
            (-4 * cos(orbit_ang_vel * t) / orbit_ang_vel - ...
                3 * orbit_ang_vel * t^2 / 2) / orbit_ang_vel;
        -3 * cos(orbit_ang_vel * t), ...
            0, ...
            sin(orbit_ang_vel * t) / orbit_ang_vel,...
            -2 * cos(orbit_ang_vel * t) / orbit_ang_vel;
        -6 * orbit_ang_vel * (t - sin(orbit_ang_vel * t) / ...
            orbit_ang_vel), ...
            0,...
            2 * cos(orbit_ang_vel * t) / orbit_ang_vel,...
            4 * sin(orbit_ang_vel * t) / orbit_ang_vel - 3 * t
    ];

    cont_input_matrix = [
        0, 0; 
        0, 0;
        1/chief_mass, 0; 
        0, 1/chief_mass
    ];

    input_matrix = (int_eat(sampling_period) - int_eat(0)) * cont_input_matrix;


end

function [state_matrix, input_matrix] = get6dCwhStateAndInputMatrices(params)
% SReachTools/getCwhLtiSystem/get4dCwhStateAndInputMatrices: Get 6-d CWH 
% matrices
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
%   state_matrix - CWH state matrix
%   input_matrix - CWH input matrix
%
% =============================================================================
%
%   This function is part of the Stochastic Optimal Control Toolbox.
%   License for the use of this function is given in
%        https://github.com/abyvinod/SReachTools/blob/master/LICENSE
% 

    % redefine important variable for simplicity
    sampling_period = params.sampling_period;
    disc_orbit_dist = params.disc_orbit_dist;
    orbit_ang_vel = params.orbit_ang_vel;
    chief_mass = params.chief_mass;

    % System matrix construction
    state_matrix = [ 
        4 - 3 * cos(disc_orbit_dist), ...
            0, ...
            0, ...
            (1/orbit_ang_vel) * sin(disc_orbit_dist), ...
            (2/orbit_ang_vel) * (1 - cos(disc_orbit_dist)), ...
            0; ...
        6 * (sin(disc_orbit_dist) - disc_orbit_dist), ...
            1, ...
            0, ...
            -(2/orbit_ang_vel) * (1 - cos(disc_orbit_dist)), ...
            (1/orbit_ang_vel) * (4*sin(disc_orbit_dist) - 3*disc_orbit_dist), ...
            0; ...
        0, ...
            0, ...
            cos(disc_orbit_dist), ...
            0, ...
            0, ...
            (1/orbit_ang_vel) * sin(disc_orbit_dist); ...
        3 * orbit_ang_vel * sin(disc_orbit_dist), ...
            0, ...
            0, ...
            cos(disc_orbit_dist), ...
            -2 * sin(disc_orbit_dist), ...
            0; ...
        -6 * orbit_ang_vel * (1 - cos(disc_orbit_dist)), ...
            0, ...
            0, ...
            -2 * sin(disc_orbit_dist), ...
            4 * cos(disc_orbit_dist) - 3, ...
            0; ...
        0, ...
            0, ...
            -orbit_ang_vel * sin(disc_orbit_dist), ...
            0, ...
            0, ...
            cos(disc_orbit_dist);
    ];

    % Input matrix construction
    cont_input_matrix = [0, 0, 0; ...
        0, 0, 0; ...
        0, 0, 0; ...
        1/chief_mass, 0, 0; ...
        0, 1/chief_mass, 0; ...
        0, 0, 1/chief_mass];

    int_eat = @(t) [
        4*t - 3 * sin(orbit_ang_vel * t)/orbit_ang_vel, ...
            0, ...
            0, ...
            -(1/orbit_ang_vel^2) * cos(orbit_ang_vel * t), ...
            (2/orbit_ang_vel) * (t - sin(orbit_ang_vel * t)/orbit_ang_vel), ...
            0; ...
        6 * (-cos(orbit_ang_vel * t)/orbit_ang_vel - orbit_ang_vel*t^2/2), ...
            t, ...
            0, ...
            -(2/orbit_ang_vel) * (t - sin(orbit_ang_vel * t)/orbit_ang_vel), ...
            (1/orbit_ang_vel) * (4*cos(orbit_ang_vel * t)/orbit_ang_vel - ...
                3*orbit_ang_vel*t^2/2), ...
            0; ...
        0, ...
            0, ...
            sin(orbit_ang_vel * t)/orbit_ang_vel, ...
            0, ...
            0, ...
            (1/orbit_ang_vel) * -cos(orbit_ang_vel * t)/orbit_ang_vel; ...
        3 * orbit_ang_vel * -cos(orbit_ang_vel * t)/orbit_ang_vel, ...
            0, ...
            0, ...
            sin(orbit_ang_vel * t)/orbit_ang_vel, ...
            2 * cos(orbit_ang_vel * t)/orbit_ang_vel, ...
            0; ...
        -6 * orbit_ang_vel * (t - sin(orbit_ang_vel * t)/orbit_ang_vel), ...
            0, ...
            0, ...
            2 * cos(orbit_ang_vel * t)/orbit_ang_vel, ...
            4 * -sin(orbit_ang_vel * t)/orbit_ang_vel - 3*t, ...
            0; ...
        0, ...
            0, ...
            orbit_ang_vel * cos(orbit_ang_vel * t)/orbit_ang_vel, ...
            0, ...
            0, ...
            sin(orbit_ang_vel * t)/orbit_ang_vel;
    ];

    input_matrix = (int_eat(sampling_period) - int_eat(0)) * cont_input_matrix;

end
