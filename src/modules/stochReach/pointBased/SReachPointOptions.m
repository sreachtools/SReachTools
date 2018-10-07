function options = SReachPointOptions(prob_str, method_str, varargin)
% Create user-specifiable options for use with SReachPoint()
% =============================================================================
%
% SReachPointOptions creates a MATLAB struct that contains user-specifiable
% options that may be used with SReachPoint
%
% =============================================================================
%
%   options = SReachPointOptions(prob_str, method_str, varargin)
% 
% Inputs:
% -------
%   prob_str    - String specifying the problem of interest. For each case, we
%                 compute the optimal value function that maps initial states
%                 to different maximal reach probabilities
%                     1. 'first' : Stay within the safety_tube and reach the
%                                  target set early if possible
%                     2. 'term' : Stay within the safety_tube
%   method_str  - Solution technique to be used (user-specifiable
%                 options associated with each technique is enumerated)
%                     'chance-open'  -- Convex chance-constrained approach for
%                                       an open-loop controller synthesis
%                                       1. pwa_accuracy: 
%                                               Accuracy of the piecewise affine 
%                                               approximation of norminvcdf
%                                               used
%                     'chance-affine'-- Convex chance-constrained approach for
%                                       an affine controller synthesis
%                                       1. verbose: Verbosity of the 
%                                               implementation (feedback for the 
%                                               user)
%                                       2. desired_accuracy: Accuracy of the
%                                               result
%                                       3. pwa_accuracy: 
%                                               Accuracy of the piecewise affine 
%                                               approximation of norminvcdf
%                                               used
%                                       4. max_input_viol_prob:
%                                               Probabilistic relaxation of the 
%                                               hard input constraints
%                                       5. bisect_lb: Bisection lower bound
%                                       6. bisect_lb: Bisection upper bound
%                                       Difference-of-convex parameters: 
%                                       7. tau_initial: Initialization of the 
%                                               slack multiplier
%                                       8. scaling_tau: Scaling factor to the 
%                                               slack multiplier
%                                       9. tau_max: Maximum value for the 
%                                               scaling factor
%                                      10. iter_max: Maximum number of
%                                               iterations for the difference of
%                                               convex iterative algorithm
%                                      11. dc_conv_tol: Tolerance for exiting 
%                                               the iterative algorithm
%                                      12. slack_tol: Tolerance within which the
%                                               the slack variable is declared 
%                                               to be equal to the norm they are 
%                                               replacing
%                     'genzps-open'  -- Genz's algorithm + Patternsearch
%                                       1. desired_accuracy: 
%                                               Accuracy of Gaussian
%                                               integral => Accuracy of the
%                                               result
%                                       2. PSoptions: 
%                                               MATLAB struct generated
%                                               using psoptimset()
%                     'scenario-open'-- Scenario-based 
%
% Outputs:
% --------
%   options     - Collection of user-specified options for 'chance-affine'
%                 (Matlab struct created using SReachPointOptions)
%
% See also SReachPoint.
%
% Notes:
% * SReachPoint() will call this function internally using the default
%     values if SReachPointOptions()-based options is not explicitly provided
%     to SReachPoint().
% ============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 
%

    valid_prob = {'first','term'};
    valid_method= {'chance-open','chance-affine','genzps-open','scenario-open'};

    % Input parsing
    inpar = inputParser();
    inpar.addRequired('prob_str', @(x) any(validatestring(x,valid_prob)));
    inpar.addRequired('method_str', @(x) any(validatestring(x,valid_method)));

    switch lower(method_str)
        case 'genzps-open'            
            inpar.addParameter('desired_accuracy',1e-3, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0}));
            inpar.addParameter('PSoptions',psoptimset('display','off'));
            
            % Ensure that patternsearch is installed
            v = ver;            
            has_patternsearch = any(strcmp(cellstr(char(v.Name)),...
                'Global Optimization Toolbox'));
            if ~has_patternsearch
                exc = SrtSetupError(['SReachPoint with ''genzps-open'' ',...
                    'option needs MATLAB''s Global Optimization Toolbox.']);
                throw(exc);
            end
        case 'chance-open'
            % Accuracy of piecewise-affine approximation of norminvcdf
            inpar.addParameter('pwa_accuracy',1e-3, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0}));
        case 'chance-affine'
            % Verbosity of the implementation
            inpar.addParameter('verbose', 0, @(x)...
                validateattributes(x, {'numeric'}, {'scalar', 'integer',...
                    '>=',0,'<=',2}));            
            % Bisection tolerance
            inpar.addParameter('desired_accuracy',1e-2, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0}));
            % Accuracy of piecewise-affine approximation of norminvcdf
            inpar.addParameter('pwa_accuracy',1e-2, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0}));
            % Probabilistic relaxation of the hard input constraints
            inpar.addParameter('max_input_viol_prob',1e-2, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0}));
            % Bisection bounds
            inpar.addParameter('bisect_lb',0, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>=',0,'<=', 1}));
            inpar.addParameter('bisect_ub',1, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>=',0,'<=', 1}));
            % Difference-of-convex: Initialization of the slack multiplier
            inpar.addParameter('tau_initial',1, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0}));
            % Difference-of-convex: Scaling factor to the slack multiplier
            inpar.addParameter('scaling_tau',2, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0}));
            % Difference-of-convex: Max scaling factor
            inpar.addParameter('tau_max',1e5, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0}));
            % Difference-of-convex: Max iterations
            inpar.addParameter('iter_max',20, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0}));
            % Difference-of-convex: Exit condition tolerance for dc iterations
            inpar.addParameter('dc_conv_tol',1e-4, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0}));
            % Difference-of-convex: When is the slack variable declared to
            % be equal to the norm they were replacing?
            inpar.addParameter('slack_tol',1e-6, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0}));
    end
    inpar.parse(prob_str, method_str, varargin{:});
    options = inpar.Results;
end
