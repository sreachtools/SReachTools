function overapprox_set = getSReachLagOverapprox(sys, target_tube, ...
    scaled_disturbance)
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
% Inputs:
% -------
%   sys          - LtiSystem object
%   target_tube  - Tube object 
%   disturbance  - Polyhedron object (bounded disturbance set)
%
% Outputs:
% --------
%   overapprox_set - Polyhedron object for the overapproximation of the 
%                    stochastic reach set
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
        inpar.parse(sys, target_tube, scaled_disturbance);
    catch cause_exc
        exc = SrtInvalidArgsError.withFunctionName();
        exc = addCause(exc, cause_exc);
        throwAsCaller(exc);
    end
    
    tube_length = length(target_tube);
    if sys.islti()
        inverted_state_matrix = inv(sys.state_mat);
    else
        throw(SrtInternalError('LtvSystem development is on going!'));

    end

    effective_target_tube = repmat(Polyhedron(), tube_length, 1);
    effective_target_tube(end) = target_tube(end);
    if tube_length > 1
        for itt = tube_length-1:-1:1
            current_time = itt - 1;
            if ~sys.islti()
                inverted_state_matrix = inv(sys.state_mat(current_time));
            end
            
            if scaled_disturbance.isEmptySet
                % No augmentation
                new_target = effective_target_tube(itt+1);
            else
                % Compute a new target set for this iteration that is robust to 
                % the disturbance
                new_target = effective_target_tube(itt+1) + (-scaled_disturbance);
            end

            % One-step backward reach set
            one_step_backward_reach_set = inverted_state_matrix * ...
                (new_target + (-sys.input_mat(current_time) * sys.input_space));

            % Guarantee staying within target_tube by intersection
            effective_target_tube(itt) = intersect(...
                one_step_backward_reach_set, target_tube(itt));
        end
    end
    overapprox_set = effective_target_tube(1);
end
