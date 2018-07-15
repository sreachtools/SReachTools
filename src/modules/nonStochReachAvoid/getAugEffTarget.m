function aug_eff_target = getAugEffTarget(sys, target_tube_with_tZero, disturbance)
% SReachTools/getAugEffTarget
% ============================================================================
%
% This function will compute the augmented effect target via the algorithm in
% the paper:
%      [[Will fill out this once paper is actually submitted]]
%
% Usage: See examples/lagrangianApproximations.m
%   
% ============================================================================
%
% Inputs:
% -------
%   sys          - LtiSystem object
%   target_tube_with_tZero  - Cell array of Polyhedron objects 
%   disturbance  - Polyhedron object (bounded disturbance set)
%
% Outputs:
% --------
%   aug_eff_target - Polyhedron object for the augmented effective
%                                target set
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
    inpar.addRequired('target_tube_with_tZero', @(x) validateattributes(x, ...
        {'TargetTube'}, {'nonempty'}));
    inpar.addRequired('disturbance', @(x) validateattributes(x, ...
        {'Polyhedron'}, {'nonempty'}));
    
    try
        inpar.parse(sys, target_tube_with_tZero, disturbance);
    catch cause_exc
        exc = SrtInvalidArgsError.withFunctionName();
        exc = addCause(exc, cause_exc);
        throwAsCaller(exc);
    end
    
    tube_length = length(target_tube_with_tZero);
    if sys.islti()
        inverted_state_matrix = inv(sys.state_mat);
    end

    effective_target_tube = repmat(Polyhedron(), tube_length, 1);
    effective_target_tube(end) = target_tube_with_tZero(end);
    if tube_length > 1
        for tube_indx = tube_length - 1:-1:1
            current_time = tube_indx - 1;
            if ~sys.islti()
                inverted_state_matrix = inv(sys.state_mat(current_time));
            end
            
            if disturbance.isEmptySet
                % No augmentation
                new_target = effective_target_tube(tube_indx+1);
            else
                % Compute a new target set for this iteration that is robust to 
                % the disturbance
                new_target = effective_target_tube(tube_indx+1) + ...
                    (-sys.dist_mat(current_time) * disturbance);
            end

            % One-step backward reach set
            one_step_backward_reach_set = inverted_state_matrix * ...
                (new_target + (-sys.input_mat(current_time) * sys.input_space));

            % Guarantee staying within target_tube_with_tZero by intersection
            effective_target_tube(tube_indx) = intersect(...
                one_step_backward_reach_set, target_tube_with_tZero(tube_indx));
        end
    end
    aug_eff_target = effective_target_tube(1);
end
