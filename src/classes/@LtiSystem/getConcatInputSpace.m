function [concatenated_input_space_A, concatenated_input_space_b] = ...
                                         getConcatInputSpace(sys,...
                                                                   time_horizon)
% SReachTools/LtiSystem/getConcatInputSpace: Get concatenated input space 
% matrices
% ============================================================================
% 
% Computes the input_space^{time_horizon} corresponding to a given set, which
% is the set of admissible open-loop control polices
%
% Usage:
% ------
%
% % Compute the (matrix form) set of admissible open-loop control policies given
% % a LtiSystem and a time horizon
%
% sys = LtiSystem(...
%     'StateMatrix', eye(2), ...
%     'InputMatrix', ones(2,1), ...
%     'InputSpace', Polyhedron('lb', -umax, 'ub', umax));
% time_horizon = 10;
% [concatenated_input_space_A, concatenated_input_space_b] = ...
%                                     getConcatInputSpace(sys,...
%                                                               time_horizon);
% 
% ============================================================================
%
% [concatenated_input_space_A, ...
%     concatenated_input_space_b] = getConcatInputSpace(sys, time_horizon)
% 
% Inputs:
% -------
%   sys                        - An object of LtiSystem class 
%   time_horizon               - Time horizon
%
% Outputs:
% --------
%   concatenated_input_space_A - State matrix for concatenated space
%   concatenated_input_space_b - Input matrix for concatenated space
%
% =============================================================================
%
% This function is part of the Stochastic Optimal Control Toolbox.
% License for the use of this function is given in
%      https://github.com/abyvinod/SReachTools/blob/master/LICENSE
% 
%

    %% Input handling
    % Ensure that the system has a non-empty input space
    assert(~sys.input_space.isEmptySet,...
           'SReachTools:invalidArgs',...
           'Expected a non-empty polyhedral input space');
    % Ensure that time horizon is a scalar and positive
    assert( isscalar(time_horizon) && time_horizon > 0,...
           'SReachTools:invalidArgs',...
           'Expected a scalar positive time_horizon');

    %% Construction of the concatenated input space
    concatenated_input_space_A = kron(eye(time_horizon), sys.input_space.A);
    concatenated_input_space_b = kron(ones(time_horizon,1), sys.input_space.b);
end
