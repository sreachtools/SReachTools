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
%                     1. 'term' : Stay within the safety_tube
%   method_str  - Solution technique to be used (user-specifiable
%                 options associated with each technique is enumerated)
%                     'chance-open'  -- Convex chance-constrained approach for
%                                       an open-loop controller synthesis
%                                       1. pwa_accuracy: Accuracy of the
%                                               piecewise affine approximation
%                                               of norminvcdf used [Default:
%                                               1e-3]
%                     'genzps-open'  -- Genz's algorithm + Patternsearch
%                                       1. desired_accuracy: Accuracy of
%                                               Gaussian integral => Accuracy of
%                                               the result [Default: 5e-2]
%                                       2. PSoptions: MATLAB struct generated
%                                               using psoptimset()
%                                       3. thresh: An upper bound on useful
%                                               reach probability. This can be 
%                                               used to specify reach_prob >= 
%                                               thresh \in (0,1] See notes 
%                                               [Default: 1]
%                     'particle-open'-- Particle control approach that uses
%                                       mixed-integer linear programming
%                                       1. n_particles: Number of particles to
%                                               use [Default: 100] | Must
%                                               be less than max_particles.
%                                       2. bigM: A large positive constant value
%                                               that is used in the mixed
%                                               integer formulation [Default:
%                                               5000]
%                                       3. verbose: Verbosity of the 
%                                               implementation (feedback for the
%                                               user) | Takes values from 0 to 2
%                                               [Default: 0]
%                                       4. max_particles: Maximum particles
%                                               permitted. This bound is
%                                               used to throw a pre-emptive
%                                               error for very demanding
%                                               problem requirements
%                                               [Default: 200]
%                     'voronoi-open' -- Voronoi-based undersampling of particle
%                                       control approach to compute open loop
%                                       1. failure_risk: Risk of the
%                                               probabilistic overapproximation
%                                               bound failing [Default: 1e-4]
%                                       2. max_overapprox_err: Maximum
%                                               overapproximation error
%                                               (probabilistically) tolerable up
%                                               to the failure_risk
%                                               [Default: 1e-2]
%                                       3. n_kmeans: Number of kmeans cluster
%                                               points/ Voronoi centers
%                                               [Default: 30]
%                                       4. bigM: A large positive constant value
%                                               used in the mixed integer 
%                                               formulation [Default: 100]
%                                       5. verbose: Verbosity of the 
%                                               implementation (feedback for the
%                                               user) | Takes values from 0 to 2
%                                               [Default: 0]
%                                       6. max_particles: Maximum particles
%                                               permitted. This bound is
%                                               used to throw a pre-emptive
%                                               error for very demanding
%                                               problem requirements
%                                               [Default: 1e5]
%                     'chance-affine'-- Difference-of-convex chance-constrained 
%                                       approach for affine controller synthesis
%                                       1. [MUST HAVE] max_input_viol_prob:
%                                               Probabilistic relaxation of the
%                                               hard input constraints 
%                                               [Default: 1e-2]
%                                       2. verbose: Verbosity of the 
%                                               implementation (feedback for the
%                                               user) | Takes values from 0 to 2
%                                               [Default: 0]
%                                       3. pwa_accuracy: Accuracy of the
%                                               piecewise affine approximation
%                                               of norminvcdf used [Default:
%                                               1e-3]
%                                       Difference-of-convex parameters: 
%                                       4. tau_initial: Initialization of the 
%                                               slack multiplier [Default: 1]
%                                       5. scaling_tau: Scaling factor to the 
%                                               slack multiplier [Default: 2]
%                                       6. tau_max: Maximum value for the 
%                                               scaling factor [Default: 1e5]
%                                       7. iter_max: Maximum number of
%                                               iterations for the difference of
%                                               convex iterative algorithm 
%                                               [Default: 200]
%                                       8. dc_conv_tol: Tolerance for exiting 
%                                               the iterative algorithm
%                                               [Default: 1e-4]
%                                       9. slack_tol: Tolerance for the sum
%                                               of slack vars for penalty DC
%                                               [Default: 1e-8]
%                     'chance-affine-uni'
%                                    -- Uniform risk allocation for affine 
%                                       controller synthesis
%                                       1. [MUST HAVE] max_input_viol_prob:
%                                               Probabilistic relaxation of the
%                                               hard input constraints 
%                                               [Default: 1e-2]
%                                       2. verbose: Verbosity of the 
%                                               implementation (feedback for the
%                                               user) | Takes values from 0 to 1
%                                               [Default: 0]
%                                       3. state_bisect_tol: Bisection
%                                               tolerance for the state 
%                                               constraint [Default: 1e-3]
%                                       4. input_bisect_tol: Bisection
%                                               tolerance for the input 
%                                               constraint [Default: 1e-3]
%
% Outputs:
% --------
%   options     - Collection of user-specified options for given method_str
%
% See also SReachPoint.
%
% Notes:
% * SReachPoint() will call this function internally using the default
%   values if SReachPointOptions()-based options is not explicitly provided
%   to SReachPoint().
% * To specify a desired set of samples V to use when undersampling in 
%   voronoi-X, set the undersampling fraction to be very small (say 1e-4/1e-5) 
%   and set min_samples to V.
% * For sampling-based approaches, we impose a heuristic maximum allowed
%   particles. This bound is used to throw a pre-emptive error for very 
%   demanding problem requirements. The user MAY MODIFY it to allow the
%   computations beyond the limit.
% * max_input_viol_prob has been left as addParameter instead of
%   addRequired, so that default values may be specified.
% * For 'genzps-open', options.thresh permits early termination if instead
%   of maximal reach probability is of interest. While this function is for 
%   maximal reach probability, this optional feature permits the reuse of
%   code in SReachSet.
% * In 'genzps-open', desired accuracy is the farthest lower bound on the
%   confidence interval acceptable. In order to remain conservative,
%   RandomVector/getProbPolyhedron subtracts desired_accuracy from the result to
%   yield an underapproximation. For higher desired_accuracy, the result may be
%   more conservative but faster. For lower desired_accuracy, the result may
%   take more time.
% ============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://sreachtools.github.io/license/
% 
%

    valid_prob = {'term'};
    valid_method= {'chance-open','chance-affine','genzps-open',...
        'particle-open','voronoi-open','chance-affine-uni'};
    
    % Input parsing
    inpar = inputParser();
    inpar.addRequired('prob_str', @(x) any(validatestring(x,valid_prob)));
    inpar.addRequired('method_str', @(x) any(validatestring(x,valid_method)));

    validatestring(method_str,valid_method);
    switch lower(method_str)
        case 'genzps-open'            
            % Ensure that patternsearch is installed
            v = ver;            
            has_patternsearch = any(strcmp(cellstr(char(v.Name)), ...
                'Global Optimization Toolbox'));
            if ~has_patternsearch
                exc = SrtSetupError(['SReachPoint with ''genzps-open'' ', ...
                    'option needs MATLAB''s Global Optimization Toolbox.']);
                throw(exc);
            end
            inpar.addParameter('desired_accuracy', 5e-2, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>=',1e-2,'<',1}));
            inpar.addParameter('PSoptions',psoptimset('display','off'));
            inpar.addParameter('thresh', 1, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0, '<=', 1}));
        case 'chance-open'
            % Accuracy of piecewise-affine approximation of norminvcdf
            inpar.addParameter('pwa_accuracy',1e-3, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0}));
        case 'particle-open'
            % Number of particles to be used for approximation
            inpar.addParameter('n_particles', 100, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0}));
            % BigM notation requires a large value
            inpar.addParameter('bigM', 100, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0}));
            % Verbosity of the implementation
            inpar.addParameter('verbose', 0, @(x)...
                validateattributes(x, {'numeric'}, {'scalar', 'integer', ...
                    '>=',0,'<=',2}));   
            % Maximum number of particles allowed for tractable computation
            inpar.addParameter('max_particles', 200, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','integer','>',0}));
         case 'voronoi-open'
            % Risk of the probabilistic overapproximation bound failing
            inpar.addParameter('failure_risk', 1e-4, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0}));
            % Maximum overapproximation error (probabilistically) tolerable
            inpar.addParameter('max_overapprox_err', 1e-2, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0}));
            % Number of kmeans cluster points/ Voronoi centers
            inpar.addParameter('n_kmeans', 30, @(x)...
                validateattributes(x, {'numeric'}, ...
                    {'scalar','>',0, 'integer'}));
            % Maximum number of particles allowed for tractable computation
            inpar.addParameter('max_particles', 1e5, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','integer','>',0}));
            % BigM notation requires a large value
            inpar.addParameter('bigM', 100, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0}));
            % Verbosity of the implementation
            inpar.addParameter('verbose', 0, @(x)...
                validateattributes(x, {'numeric'}, {'scalar', 'integer', ...
                    '>=',0,'<=',2}));                        
        case 'chance-affine'
            % Probabilistic relaxation of the hard input constraints
            inpar.addParameter('max_input_viol_prob',1e-2, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0,'<',1}));
            % Verbosity of the implementation
            inpar.addParameter('verbose', 0, @(x)...
                validateattributes(x, {'numeric'}, {'scalar', 'integer', ...
                    '>=',0,'<=',2}));            
            % Accuracy of piecewise-affine approximation of norminvcdf
            inpar.addParameter('pwa_accuracy',1e-3, @(x)...
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
        case 'chance-affine-uni'
            % Probabilistic relaxation of the hard input constraints
            inpar.addParameter('max_input_viol_prob',1e-2, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0,'<',1}));
            % Verbosity of the implementation
            inpar.addParameter('verbose', 0, @(x)...
                validateattributes(x, {'numeric'}, {'scalar', 'integer', ...
                    '>=',0,'<=',1}));            
            % Bisection tolerance for the state chance constraints
            inpar.addParameter('state_bisect_tol', 1e-3, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0,'<',1}));
            % Bisection tolerance for the input chance constraints
            inpar.addParameter('input_bisect_tol', 1e-3, @(x)...
                validateattributes(x, {'numeric'}, {'scalar','>',0,'<',1}));           
    end
    inpar.parse(prob_str, method_str, varargin{:});
    options = inpar.Results;
 
    % Check for must-have parameters
    switch lower(method_str)
        case {'chance-affine','chance-affine-uni'}
            if any(strcmpi(inpar.UsingDefaults, 'max_input_viol_prob'))
                 throwAsCaller(SrtInvalidArgsError(['Expected ', ...
                     'max_input_viol_prob, the maximum allowed likelihood ', ...
                     'of violating the input constraints.']));
            end
        case 'particle-open' 
            if options.n_particles > options.max_particles
                throwAsCaller(SrtInvalidArgsError(...
                     sprintf(['Particles required (%d) > maximum allowed ',...
                        'particles (%d).'], options.max_particles,...
                        options.n_particles)));
            end
    end
end
