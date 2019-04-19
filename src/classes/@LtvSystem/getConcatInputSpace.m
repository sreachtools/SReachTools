function [concat_input_space_A, concat_input_space_b] = ...
                                         getConcatInputSpace(sys, ...
                                                             time_horizon)
% Get half space representation of the concatenated (polytopic) input space 
% for the given time horizon
% ============================================================================
% 
% Computes the input_space^{time_horizon} corresponding to a given set, which
% is the set of admissible open-loop control polices. This function computes the
% half-space representation of the cartesian products of polytopic input spaces.
%
% Usage:
% ------
%
% % Compute the (matrix form) set of admissible open-loop control policies given
% % a LtvSystem and a time horizon
%
% sys = LtvSystem(...
%     'StateMatrix', eye(2), ...
%     'InputMatrix', ones(2,1), ...
%     'InputSpace', Polyhedron('lb', -umax, 'ub', umax));
% time_horizon = 10;
% [concat_input_space_A, concat_input_space_b] = ...
%                                             getConcatInputSpace(sys, ...
%                                                                 time_horizon);
% 
% ============================================================================
%
% [concat_input_space_A, concat_input_space_b] =...
%                                              getConcatInputSpace(sys, ...
%                                                                  time_horizon)
% 
% Inputs:
% -------
%   sys                  - An object of LtvSystem class 
%   time_horizon         - Time horizon
%
% Outputs:
% --------
%   concat_input_space_A, concat_input_space_b 
%                        - Concatenated input space (Halfspace representation)
%
% =============================================================================
%
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://sreachtools.github.io/license/
% 
%

    %% Input handling
    % Ensure that the system has a non-empty input space
    if sys.input_space.isEmptySet
        throwAsCaller(SrtInvalidArgsError(['Expected a non-empty ', ...
            'polyhedral input space']));
    end

    % Ensure that time horizon is a scalar and positive
    if ~isscalar(time_horizon) || time_horizon <= 0
        throwAsCaller(SrtInvalidArgsError(['Expected a scalar ', ...
            'positive time_horizon']));
    end

    %% Construction of the concatenated input space (input_space^{time_horizon})
    concat_input_space_A = kron(eye(time_horizon), sys.input_space.A);
    concat_input_space_b = kron(ones(time_horizon,1), sys.input_space.b);
end
