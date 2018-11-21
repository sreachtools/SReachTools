function varargout = getSReachLagOverapprox(sys, target_tube, disturbance_set)
% Get the overapproximation of the stoch reach set
% ============================================================================
%
% This function will compute the overapproximation of the stochastic reach
% set via Algorithm 2 in
% 
%      J. D. Gleason, A. P. Vinod, and M. M. K. Oishi. 2018. Lagrangian 
%      Approximations for Stochastic Reachability of a Target Tube. 
%      online. (2018). https://arxiv.org/abs/1810.07118
%
% Usage: See examples/lagrangianApproximations.m
%   
% ============================================================================
%
% [overapprox_set, overapprox_tube] = getSReachLagUnderapprox(sys,...
%       target_tube, disturbance_set)
%
% Inputs:
% -------
%   sys             - LtiSystem object
%   target_tube     - Tube object 
%   disturbance_set - Polyhedron object (bounded disturbance set)
%
% Outputs:
% --------
%   overapprox_set - Polyhedron object for the overapproximation of the 
%                    stochastic reach set
%   overapprox_tube- [Optional] Tube comprising of an overapproximation of the
%                    stochastic reach sets across the time horizon
%
% Notes:
% * From computational geometry, intersections and Minkowski differences are
%   best performed in facet representation and Minkowski sums are best
%   performed in vertex representation. However, since in this computation,
%   all three operations are required, scalability of the algorithm is severly
%   hampered, despite theoretical elgance.
%
% ============================================================================
%
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
%

    inpar = inputParser();
    inpar.addRequired('sys', @(x) validateattributes(x, ...
        {'LtiSystem', 'LtvSystem'}, {'nonempty'}));
    inpar.addRequired('target_tube', @(x) validateattributes(x, ...
        {'Tube'}, {'nonempty'}));
    inpar.addRequired('disturbance', @(x) validateattributes(x, ...
        {'Polyhedron'}, {'nonempty'}));
    
    try
        inpar.parse(sys, target_tube, disturbance_set);
    catch cause_exc
        exc = SrtInvalidArgsError.withFunctionName();
        exc = addCause(exc, cause_exc);
        throwAsCaller(exc);
    end
    
    tube_length = length(target_tube);
    if sys.islti()
        inverted_state_matrix = inv(sys.state_mat);
        minus_bu = (-sys.input_mat) * sys.input_space;
        minus_scaled_dist_set = (-sys.dist_mat) * disturbance_set;
    end

    effective_target_tube = repmat(Polyhedron(), tube_length, 1);
    effective_target_tube(end) = target_tube(end);
    if tube_length > 1
        for itt = tube_length-1:-1:1
            current_time = itt - 1;
            if sys.isltv()
                % Overwrite the following parameters with their
                % time-varying counterparts
                inverted_state_matrix = inv(sys.state_mat(current_time));
                minus_bu = (-sys.input_mat(current_time)) * sys.input_space;
                minus_scaled_dist_set = (-sys.dist_mat(current_time)) *...
                    disturbance_set;
            end
            
            if disturbance_set.isEmptySet
                % No augmentation
                new_target = effective_target_tube(itt+1);
            else
                % Compute a new target set for this iteration that is robust to 
                % the disturbance
                new_target = effective_target_tube(itt+1)+minus_scaled_dist_set;
            end

            % One-step backward reach set
            one_step_backward_reach_set = inverted_state_matrix * ...
                (new_target + minus_bu);

            % Guarantee staying within target_tube by intersection
            effective_target_tube(itt) = intersect(...
                one_step_backward_reach_set, target_tube(itt));
        end
    end
    varargout{1} = effective_target_tube(1);
    if nargout > 1
        varargout{2} = effective_target_tube;
    end
end
