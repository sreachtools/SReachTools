function [sys, heading_vec] = getDubinsCarLtv(type,...
    turning_rate_seq,...
    initial_heading,...
    sampling_time,...
    varargin)
% SReachTools/systems/getDubinsCarLtv: 
% ============================================================================
% 
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
%        https://github.com/abyvinod/SReachTools/blob/master/LICENSE
% 
%

    heading_vec = initial_heading + sampling_time * cumsum([0;turning_rate_seq]);

    switch(lower(type))
        case 'add-dist'
            % TODO: Check if it is a proper Gaussian variable            
            assert(isa(varargin{1}, 'Polyhedron') && varargin{1}.Dim == 1, ...
               'SReachTools:invalidArgs',...
               'Input space is a polytope');           
           
            sys = LtvSystem('StateMatrix', @(t) eye(2), ...
                'InputMatrix', @(t) sampling_time * [cos(heading_vec(t+1)); ...
                                                     sin(heading_vec(t+1))], ...
                'InputSpace', varargin{1},...
                'DisturbanceMatrix', varargin{2},...
                'Disturbance', varargin{3});            
        case 'vel-dist'
            % TODO: Throw error if more than 1 varargin
            % TODO: Check if it is a proper Gaussian variable
            sys = LtvSystem('StateMatrix', @(t) eye(2), ...
                'DisturbanceMatrix', @(t) sampling_time * [cos(heading_vec(t+1)); ...
                                                           sin(heading_vec(t+1))], ...
                'Disturbance', varargin{1});
        otherwise
    end    
end