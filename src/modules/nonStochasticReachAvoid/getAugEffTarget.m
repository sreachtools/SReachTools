function aug_eff_target = getAugEffTarget(sys, ...
                                          target_tube, ...
                                          disturbance)
% SReach/getAugEffTarget
% ============================================================================
%
% This function will compute the augmented effect target via the algorithm in
% the paper:
%      [[Will fill out this once paper is actually submitted]]
%
% Usage
% -----
% Will create example to demonstrate use
%   
% ============================================================================
%
% Inputs:
% -------
%   sys          - LtiSystem object
%   target_tube  - Cell array of Polyhedron objects 
%   disturbance  - Polyhedron object (bounded disturbance set)
%
% Outputs:
% --------
%   aug_eff_target - Polyhedron object for the augmented effective
%                                target set
%
% Notes:
%   - From computational geometry, intersections and Minkowski differences are
%     best performed in facet representation and Minkowski sums are best
%     performed in vertex representation. However, since in this computation,
%     all three operations are required, scalability of the algorithm is severly
%     hampered, despite theoretical elgance.
%
% ============================================================================
%
%   This function is part of the Stochastic Optimal Control Toolbox.
%   License for the use of this function is given in
%        https://github.com/abyvinod/SReach/blob/master/LICENSE
%

    % validate the inputs
    validateattributes(sys, {'LtiSystem'}, {'nonempty'});
    validateattributes(target_tube, {'cell'}, {'nonempty'});
    validateattributes(disturbance, {'Polyhedron'}, {'nonempty'});
    
    % validate that all elements of the target_tube are polyhedron
    for i = 1:length(target_tube)
        validateattributes(target_tube{i}, {'Polyhedron'}, {'nonempty'});
    end
    
    horizon_length = length(target_tube);
    inverted_state_matrix = inv(sys.state_matrix);
    eff_target_temp = target_tube{horizon_length};
    if horizon_length > 1
        for i = 1:horizon_length - 1
            if disturbance.isEmptySet
                % No requirement of robustness
                new_target = eff_target_temp;
            else
                % Compute a new target set for this iteration that is robust to 
                % the disturbance
                new_target = eff_target_temp + ...
                             (-sys.disturbance_matrix * disturbance);
            end

            % One-step backward reach set
            one_step_backward_reach_set = inverted_state_matrix * ...
                (new_target + (-sys.input_matrix * sys.input_space));


            % Guarantee staying within target_tube by intersection
            eff_target_temp = intersect(one_step_backward_reach_set, ...
                                              target_tube{horizon_length-i});
        end
    end
    aug_eff_target = eff_target_temp;
end