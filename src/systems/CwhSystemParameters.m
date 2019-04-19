classdef CwhSystemParameters
% A MATLAB class to store/retrive the default parameters used in CWH dynamics
% with modifications, if any
% =============================================================================
%
% Notes:
% ------
% * This code and the parameters were obtained from Lesser's repeatability code
%   for the 2013 CDC paper.
% * The default parameters (with their Names to specify changes) are:
%       sampling_period : sampling period        = 20 s
%       orbital_radius  : orbital radius         = 850 + 6378.1 m
%       grav_constant   : gravitational constant = 6.673e-11
%       celes_mass      : celestial body mass    = 5.9472e24 kg
%       chief_mass      : chief mass             = 300 kg
% * Along with these parameters, the class provides these parameters that are
%   computed using the above parameters
%       grav_body       : gravitational body           = grav_constant *
%                                                           celes_mass / 1e6
%       orbit_ang_vel   : orbital angular velocity     = sqrt(grav_body /
%                                                              orbital_radius^3)
%       disc_orbit_dist : discretized orbital distance = orbit_ang_vel *
%                                                           sampling_period rad
%
% =============================================================================
%
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://sreachtools.github.io/license/
% 

    properties
        % sampling period in sec
        % Default value: 20
        sampling_period

        % Orbital radius in m
        % Default value: 7228.1
        orbital_radius

        % Universal gravitational constant  
        % Default value: 6.673e-11
        grav_constant

        % Mass of the celestial body (default is Earth)
        % Default value: 5.9472e24
        celes_mass

        % Mass of the chief kg
        % Default value: 300
        chief_mass

        % Gravitation constant for the pull of the celestial body 
        % (default Earth)
        % Set via the equation
        %   grav_body = grav_constant * celes_mass / 1000^3;
        grav_body

        % Angular velocity in the orbit
        % Set via the equatoin
        %   orbit_ang_vel = sqrt(grav_body / orbital_radius^3);
        orbit_ang_vel

        % Discretized orbital distance
        % Set via the equation
        %    disc_orbit_dist = orbit_ang_vel * sampling_period;   
        disc_orbit_dist

    end

    methods
        function obj = CwhSystemParameters(varargin)

            % input handling
            inpar = inputParser();
            inpar.PartialMatching = false;
            inpar.KeepUnmatched = true;
            inpar.addParameter('SamplingPeriod', 20, ...
                @(x) validateattributes(x, {'numeric'}, {'scalar', '>', 0}));
            inpar.addParameter('OrbitalRadius', 850+6378.1, ...
                @(x) validateattributes(x, {'numeric'}, {'scalar', '>', 0}));
            inpar.addParameter('GravConstant', 6.673e-11, ...
                @(x) validateattributes(x, {'numeric'}, {'scalar', '>', 0}));
            inpar.addParameter('CelestialMass', 5.9472e24, ...
                @(x) validateattributes(x, {'numeric'}, {'scalar', '>', 0}));
            inpar.addParameter('ChiefMass', 300, ...
                @(x) validateattributes(x, {'numeric'}, {'scalar', '>', 0}));

            try
                inpar.parse(varargin{:});
            catch err
                exc = SrtInvalidArgsError.withFunctionName();
                exc = exc.addCause(err);
                throwAsCaller(exc);
            end
            if ~isempty(fieldnames(inpar.Unmatched))
                err_mesg = 'Invalid argument(s)---';
                invalid_args = fieldnames(inpar.Unmatched);
                for indx = 1:length(invalid_args)
                    err_mesg = [err_mesg,' ',char(invalid_args{indx}),','];
                end
                err_mesg=[err_mesg, ...
                    '\b---given to CwhSystemParameters. Expected'];
                for valid_arg = inpar.Parameters
                    err_mesg = [err_mesg,' ',char(valid_arg),','];
                end
                err_mesg = strcat(err_mesg,'\b as the (optional) arguments.');
                exc = SrtInvalidArgsError(err_mesg);
                throwAsCaller(exc);
            end

            obj.sampling_period = inpar.Results.SamplingPeriod;
            obj.orbital_radius  = inpar.Results.OrbitalRadius;
            obj.grav_constant   = inpar.Results.GravConstant;
            obj.celes_mass      = inpar.Results.CelestialMass;
            obj.chief_mass      = inpar.Results.ChiefMass;

            % Gravitation constant for the pull of the celestial body 
            % (default Earth)
            obj.grav_body = obj.grav_constant * obj.celes_mass / 1000^3;

            % Angular velocity in the orbit
            obj.orbit_ang_vel = sqrt(obj.grav_body / obj.orbital_radius^3);

            % Discretized orbital distance
            obj.disc_orbit_dist = obj.orbit_ang_vel * obj.sampling_period;   

        end
    end

end
