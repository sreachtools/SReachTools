function [sys, heading_vec] = getDubinsCarLtv(type, turning_rate_seq, ...
    initial_heading, sampling_time, varargin)
% Get a LtvSystem object for the Dubins car model with known turning rate
% sequence
% ============================================================================
% 
% Usage:
% ------
% % Known turning rate sequence
% sampling_time = 0.1;                        % Sampling time
% time_horizon = 50;                          % Max number of time steps in
%                                             % simulation
% init_heading = pi/10;                       % Initial heading 
% % Create a constant turning rate sequence
% omega = pi/time_horizon/sampling_time;      
% turning_rate = omega*ones(time_horizon,1);   
% % Input space definition
% umax = 6;
% input_space = Polyhedron('lb',0,'ub',umax);
% % Disturbance matrix and random vector definition
% dist_matrix = eye(2);
% eta_dist = RandomVector('Gaussian',zeros(2,1), 0.001 * eye(2));
% 
% [sys, heading_vec] = getDubinsCarLtv('add-dist', turning_rate, ...
%   init_heading, sampling_time, input_space, dist_matrix, eta_dist);
%
% ============================================================================
% 
% sys = getDubinsCarLtv(type, turning_rate_seq, initial_heading, sampling_time,
%   varargin)
%
% Inputs:
% -------
%   type        - 
%   turning_rate_seq    - Known turning rate sequence (column vector of length
%                         time_horizon) 
%   initial_heading     - Initial heading angle
%   sampling_time       - Sampling time for the system
%   Required additional arguments for different types:
%   type: 'add-dist' (Additive disturbance with velocity as input)
%       velocity_input  - Bounds on the velocity (1-dimensional Polyhedron
%                         object)
%       dist_matrix     - Disturbance matrix for the additive disturbance 
%       dist            - Disturbance
%   type = 'vel-dist' (Velocity is the disturbance) 
%       velocity_dist   - Velocity disturbance
% 
% Outputs:
% --------
%   sys                 - LtvSystem object describing the Dubin's car
% 
% =============================================================================
%
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/abyvinod/SReachTools/blob/master/LICENSE
% 
%

    valid_types = {'add-dist','vel-dist'};
    validatestring(type, valid_types, 'getDubinsCarLtv', 'type');
    validateattributes(initial_heading, {'numeric'}, {'scalar'}, ...
        'getDubinsCarLtv', 'initial_heading');
    validateattributes(sampling_time, {'numeric'}, {'scalar','>',0}, ...
        'getDubinsCarLtv', 'sampling_time');
    validateattributes(turning_rate_seq, {'numeric'}, {'column'}, ...
        'getDubinsCarLtv', 'turning_rate_seq');

    heading_vec = initial_heading + sampling_time * ...
        cumsum([0;turning_rate_seq]);

    time_varying_matrix = @(t) sampling_time * [cos(heading_vec(t+1)); ...
                                                sin(heading_vec(t+1))];
    switch(lower(type))
        case 'add-dist'
            if length(varargin) < 3
                throwAsCaller(SrtInvalidArgsError('Too few arguments'));
            end
            velocity_input = varargin{1};
            validateattributes(velocity_input, {'Polyhedron'},...
                {'nonempty'}, 'getDubinsCarLtv',...
                'velocity_input with add-dist option');
            if velocity_input.Dim ~= 1
                throwAsCaller(SrtInvalidArgsError(['Velocity bounds ',...
                    '(velocity_input) must be one-dimensional polyhedron']));
            end
            dist_matrix = varargin{2};
            validateattributes(dist_matrix, {'numeric'},...
                {'nonempty'}, 'getDubinsCarLtv',...
                'dist_matrix with add-dist option');
            dist = varargin{3};
            validateattributes(dist, {'RandomVector', 'Polyhedron'},...
                {'nonempty'}, 'getDubinsCarLtv',...
                'dist with add-dist option');
            if length(varargin) > 3
                throwAsCaller(SrtInvalidArgsError('Too many input arguments'));
            end
            sys = LtvSystem('StateMatrix', @(t) eye(2), ...
                'InputMatrix', time_varying_matrix, ...
                'InputSpace', velocity_input, 'Disturbance', dist, ...
                'DisturbanceMatrix', dist_matrix);
        case 'vel-dist'
            if isempty(varargin)
                throwAsCaller(SrtInvalidArgsError('Too few arguments'));
            end
            velocity_dist = varargin{1};
            validateattributes(velocity_dist, {'RandomVector', 'Polyhedron'},...
                {'nonempty'}, 'getDubinsCarLtv',...
                'velocity_dist with vel-dist option');
            if length(varargin) > 1
                throwAsCaller(SrtInvalidArgsError('Too many input arguments'));
            end
            sys = LtvSystem('StateMatrix', @(t) eye(2), ...
                'DisturbanceMatrix', time_varying_matrix, ...
                'Disturbance', velocity_dist);
        otherwise
    end    
end
