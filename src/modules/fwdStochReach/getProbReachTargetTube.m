function prob = getProbReachTargetTube(sys, ...
                                       initial_state, ...
                                       target_tube, ...
                                       desired_accuracy, ...
                                       varargin)
% SReachTools/forwardStochasticReach/getProbReachTargetTube: Compute the
% probability that the state will lie in a target tube. The starting point may
% be a vector or a RandomVector object
% ============================================================================
%
% This function uses getHmatMeanCovForXSansInput to compute the forward
% stochastic reach probability density (FSRPD) of the concatenated state vector.
% Next, it evaluates the integral of the resulting Gaussian over the
% user-specified target tube (a MPT Polyhedron) using iteratedQscmvnv.
%
% Usage: See examples/forwardStochasticReachCWH.mlx
%
% ============================================================================
% 
% prob = getProbReachTargetTube(sys, ...
%                               initial_state, ...
%                               target_tube, ...
%                               desired_accuracy, ...
%                               varargin)
%
% Inputs:
% -------
%   sys              - An object of LtiSystem class 
%   initial_state    - Initial state can be a deterministic n-dimensional vector
%                      or a RandomVector object
%   target_tube      - Target tube to stay within [TargetTube object]
%   desired_accuracy - Accuracy of the integral evaluation 
%                      [Default 1e-3 otherwise]
%   input_policy     - (Required only for controlled systems) Input policy 
%
% Outputs:
% --------
%   prob             - Probability that the system accomplishes the reach-avoid
%                      objective
%
% See also iteratedQscmvnv, getFSRPDMeanCov.
%
% Notes:
% ------
% * In case, the target set is a hyper-cuboid and the state_dim < 25,
%   then use mvncdf instead.
% * The safe set and the target set must be Polyhedron objects.
% ============================================================================
%
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
%
%

    % Create concatenated target tube 
    [concat_target_tube_A, concat_target_tube_b] = target_tube.concat();
    time_horizon = length(target_tube);
    
    %% INPUT HANDLING
    % Scalar desired_accuracy
    assert(isscalar(desired_accuracy), ...
           'SReachTools:invalidArgs', ...
           'Expected a scalar value for desired_accuracy');

    % Obtain the input policy
    if sys.input_dim > 0
        assert(length(varargin) == 1, ...
               'SReachTools:invalidArgs', ...
               'Expected an input policy only for a controlled system');
        input_policy = varargin{1};
    else
        input_policy = zeros(0, 1);
    end

    % Compute H (zeros(sys.state_dim*time_horizon,1)), mean_X_sans_input,
    % cov_X_sans_input for the safety_cost_function definition
    % GUARANTEES: Gaussian-perturbed LTI system (sys) and well-defined
    % initial_state and time_horizon
    [H, mean_X_sans_input, cov_X_sans_input] = ...
        getHmatMeanCovForXSansInput(sys, initial_state, time_horizon);
        
    % Compute the reach-avoid probability for the uncontrolled system
    prob = computeReachAvoidProb( ...
                input_policy, ...
                mean_X_sans_input, ...
                cov_X_sans_input, ...
                H, ...
                concat_target_tube_A, ...
                concat_target_tube_b, ...
                desired_accuracy);
end
