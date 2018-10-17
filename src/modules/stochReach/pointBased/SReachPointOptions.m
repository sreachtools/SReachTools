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
%                                       2. pwa_accuracy: 
%                                               Accuracy of the piecewise affine 
%                                               approximation of norminvcdf
%                                               used
%                                       3. max_input_viol_prob:
%                                               Probabilistic relaxation of the 
%                                               hard input constraints
%                                       Difference-of-convex parameters: 
%                                       4. tau_initial: Initialization of the 
%                                               slack multiplier
%                                       5. scaling_tau: Scaling factor to the 
%                                               slack multiplier
%                                       6. tau_max: Maximum value for the 
%                                               scaling factor
%                                       7. iter_max: Maximum number of
%                                               iterations for the difference of
%                                               convex iterative algorithm
%                                       8. dc_conv_tol: Tolerance for exiting 
%                                               the iterative algorithm
%                                       9. slack_tol: Tolerance for the sum
%                                               of slack vars for penalty DC
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
    valid_method= {'chance-open','chance-affine','genzps-open','particle-open'};

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
        case 'particle-open'
            % Number of particles to be used for approximation
            inpar.addParameter('num_particles', 1e2, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0}));
            % BigM notation requires a large value
            inpar.addParameter('bigM', 5e3, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0}));
            % Verbosity of the implementation
            inpar.addParameter('verbose', 0, @(x)...
                validateattributes(x, {'numeric'}, {'scalar', 'integer',...
                    '>=',0,'<=',2}));                        
        case 'chance-affine'
            % Probabilistic relaxation of the hard input constraints
            inpar.addParameter('max_input_viol_prob',1e-2, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0,'<',1}));
            % Verbosity of the implementation
            inpar.addParameter('verbose', 0, @(x)...
                validateattributes(x, {'numeric'}, {'scalar', 'integer',...
                    '>=',0,'<=',2}));            
            % Accuracy of piecewise-affine approximation of norminvcdf
            inpar.addParameter('pwa_accuracy',1e-2, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0}));
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
            inpar.addParameter('iter_max',200, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0}));
            % Difference-of-convex: Exit condition tolerance for dc iterations
            inpar.addParameter('dc_conv_tol',1e-4, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0}));
            % Difference-of-convex: Slack tolerance requirements
            inpar.addParameter('slack_tol',1e-8, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0}));
    end
    inpar.parse(prob_str, method_str, varargin{:});
    options = inpar.Results;
    
    switch lower(method_str)
        case 'chance-affine'
            if any(strcmp(inpar.UsingDefaults, 'max_input_viol_prob'))
                 throwAsCaller(SrtInvalidArgsError(['Expected ',...
                     'max_input_viol_prob, the maximum allowed likelihood ',...
                     'of violating the input constraints.']));
            end
    end
end

