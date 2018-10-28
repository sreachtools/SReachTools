function options = SReachSetOptions(prob_str, method_str, varargin)
% Create user-specifiable options for use with SReachSet()
% =============================================================================
%
% SReachSetOptions creates a MATLAB struct that contains user-specifiable
% options that may be used with SReachSet
%
% =============================================================================
%
%   options = SReachSetOptions(prob_str, method_str, varargin)
% 
% Inputs:
% -------
%   prob_str    - String specifying the problem of interest. For each case, we
%                 compute the optimal value function that maps initial states
%                 to different maximal reach probabilities
%                     1. 'term' : Stay within the safety_tube
%   method_str  - Solution technique to be used; available techniques:
%                     'chance-opn', 'genzps-open', 'lag-under', 'lag-over'
%   varargin    - Additional required options for each technique, specified as
%                 Name-Value pairs. The additional required options are 
%                 specified below.
% 
%       'chance-open' : Convex chance-constrained approach for an open-loop 
%                       controller synthesis
%           1. set_of_dir_vecs
%               - Set of direction vectors shot outwards from the initial state 
%                 with maximum reach probability to identify the vertices of the
%                 underapproximative polytope for the stochastic reach set
%           2. init_safe_set_affine
%               - Affine constraints (if any) on the initial state. Must include 
%                 a translate of the affine hull of the set_of_dir_vecs | On 
%                 intersection with the safe set, it should result in a 2-D set
%           3. verbose
%               - Verbosity of the implementation {0,1}
%                  0 - No output 
%                  1 - Outputs the direction vector being analyzed and the 
%                      method used to obtain xmax under study (maximizing the 
%                      reach probability or the Chebyshev centering)
%           4. pwa_accuracy
%               - Accuracy of the piecewise affine overapproximation of the 
%                 inverse of the standard normal cumulative density function
% 
%       'genzps-open' : Genz's algorithm + Patternsearch
%            1. set_of_dir_vecs
%                - Set of direction vectors shot outwards from the initial state 
%                  with maximum reach probability to identify the vertices of 
%                  the underapproximative polytope for the stochastic reach set
%            2. init_safe_set_affine
%                - Affine constraints (if any) on the initial state. Must 
%                  include a translate of the affine hull of the 
%                  set_of_dir_vecs | On intersection with the safe set, it 
%                  should result in a 2-D set
%            3. verbose
%                - Verbosity of the implementation {0,1}
%                   0 - No output 
%                   1 - Outputs the direction vector being analyzed, the 
%                       summarized progress made by the bisection used for 
%                       line search
%            4. tol_bisect
%                - Tolerance for the bisection to terminate the line search 
%                  along a direction vector for the vertex of the 
%                  underapproximative polytope [Default 1e-2]
%            5. desired_accuracy
%                - Accuracy expected for the integral of the Gaussian random 
%                  vector X over the safety tube => Accuracy of the result 
%                  [Default 1e-3]
%            6. PSoptions
%                - MATLAB struct from psoptimset(), options for MATLAB's 
%                  patternsearch [Default psoptimset('Display', 'off')]
% 
%       'lag-over'/'lag-under' : Lagrangian-based over- and underapproximation
%            1. bound_set_method
%                - Method for obtaining the bounded set for over or 
%                  underapproximation. The available methods are: 'random', 
%                  'box', and 'load', each of which requires additional 
%                  arguments.
% 
%                a. random 
%                    - Get an approximation of the ellipsoid using random 
%                      direction choices; only usable for Gaussian-type 
%                      disturbances;
%                      **This method requires the additional following option 
%                        specifications:**
%                    i. num_dirs 
%                        - Number of directions to sample 
%                          the ellipsoid for a polytopic representation
%               b. box 
%                   - Get an n-dimensional rectangle centered at the disturbance 
%                     mean that satisfies the probability threshold.
%                      **This method requires the additional following option 
%                        specifications:**
%                   i. err_thresh - Tolerance for the bisection 
%                      algorithm that identifies the length of the 
%                      box
%               c. load 
%                   - Load a predefined polyhedron bounding set; primarily used 
%                     for comparison and repeatability testing.
%                      **This method requires the additional following option 
%                        specifications:**
%                   i. load_str - Path to the file to load. All 
%                      other inputs are IRRELEVANT for this option.
%                      Mat files to be loaded must have only contain
%                      the bounded set (Polyhedron object).
%
% Outputs:
% --------
%   options     - Collection of user-specified options for 'chance-affine'
%                 (Matlab struct created using SReachSetOptions)
%
% See also SReachSet.
%
% Notes:
% * SReachSetOptions() requires init_safe_set_affine and set_of_dir_vecs for the
%   methods 'chance-open' and 'genzps-open'.
% 
% ============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 
%

    valid_prob = {'term'};
    valid_method= {'chance-open','genzps-open','lag-under','lag-over'};
    valid_bound_method = {'load','random','box','ellipsoid'};
    % Input parsing
    inpar = inputParser();
    inpar.addRequired('prob_str', @(x) any(validatestring(x,valid_prob)));
    inpar.addRequired('method_str', @(x) any(validatestring(x,valid_method)));

    if contains(method_str,'lag-')
        inpar.addParameter('bound_set_method',[], @(x) any(validatestring(x, ...
            valid_bound_method)));
        %% Optional arguments that are made REQUIRED based on bound_set_method
        % Get number of directions for random option
        inpar.addParameter('num_dirs', 10, @(x) validateattributes(x, ...
            {'numeric'}, {'scalar', 'integer'}));
        % Get bisection error threshold for box option
        inpar.addParameter('err_thresh', 1e-3, @(x) validateattributes(x, ...
            {'numeric'}, {'scalar', 'positive'}));
        % Get load string for load option
        inpar.addParameter('load_str', ' ', @(x) validateattributes(x, ...
            {'char'}, {'nonempty'}));
    else
        % Hyperplane intersecting the stochastic reach set underapprox.
        inpar.addParameter('init_safe_set_affine', Polyhedron(), ...
            @(x) validateattributes(x, {'Polyhedron'}, {'nonempty'}));
        % Direction vectors to be used for line-search for vertices
        inpar.addParameter('set_of_dir_vecs', [], @(x) validateattributes(x, ...
            {'numeric'}, {'nonempty'}));
        % Verbosity
        inpar.addParameter('verbose',0, @(x) validateattributes(x, ...
            {'numeric'}, {'scalar','>=',0,'<=',1}));

        switch lower(method_str)
            case 'genzps-open'  
                % Ensure that patternsearch is installed
                v = ver;            
                has_fmincon = any(strcmp(cellstr(char(v.Name)), ...
                    'Global Optimization Toolbox'));
                if ~has_fmincon
                    exc = SrtSetupError(['SReachSet with ''genzps-open'' ', ...
                        'option needs MATLAB''s Global Optimization Toolbox.']);
                    throw(exc);
                end
                % Tolerance for bisection for the line search
                inpar.addParameter('tol_bisect',1e-2, @(x)...
                    validateattributes(x, {'numeric'}, {'scalar','>',0}));
                % Accuracy for Genz's algorithm to compute integral of Gaussian
                inpar.addParameter('desired_accuracy',1e-3, @(x)...
                    validateattributes(x, {'numeric'}, {'scalar','>',0}));
                % Patternsearch options
                inpar.addParameter('PSoptions',psoptimset('display','off'));
            case 'chance-open'
                % Accuracy of piecewise-affine approximation of norminvcdf
                inpar.addParameter('pwa_accuracy',1e-3, @(x)...
                    validateattributes(x, {'numeric'}, {'scalar','>',0}));
        end
    end
    
    %% Construct the options MATLAB struct
    inpar.parse(prob_str, method_str, varargin{:});
    options = inpar.Results;
    
    %% Ensure that user provided an option for the must use case
    if contains(method_str,'lag-') 
        if any(contains(inpar.UsingDefaults,'bound_set_method'))
            throw(SrtInvalidArgsError(['bound_set_method is a required ', ...
                'input for SReachSet when using ''lag-over''/''lag-under''.']));
        end
        switch(lower(options.bound_set_method))
            case 'box'
                if any(contains(inpar.UsingDefaults,'err_thresh'))
                    throw(SrtInvalidArgsError(['err_thresh (threshold for ', ...
                        'bisection) is a required input for ', ...
                        'bound_set_method: box']));
                end
            case 'random'
                if any(contains(inpar.UsingDefaults,'num_dirs'))
                    throw(SrtInvalidArgsError(['num_dirs (no. of ', ...
                        'direction vectors to create bounded polytope) is ', ...
                        'a required input for bound_set_method: random']));
                end
            case 'load'
                if any(contains(inpar.UsingDefaults,'load_str'))
                    throw(SrtInvalidArgsError(['Expected a path to ', ...
                        'matfile with bounded set (only) as input for ', ...
                        'bound_set_method: load']));
                end
        end
    elseif strcmpi(method_str,'chance-open') || ...
           strcmpi(method_str,'genzps-open')
           
        % Must have for non-lag options: init_safe_set_affine, set_of_dir_vecs
        if any(strcmp(inpar.UsingDefaults, 'init_safe_set_affine')) || ...
                any(strcmp(inpar.UsingDefaults, 'set_of_dir_vecs'))
             throw(SrtInvalidArgsError(['Expected init_safe_set_affine ', ...
                'and set_of_dir_vecs']));
        end            
    end
end
