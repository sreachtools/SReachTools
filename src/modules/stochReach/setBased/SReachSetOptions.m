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
%           1. set_of_dir_vecs      - [MUST HAVE] Set of direction vectors shot
%                                     outwards from the initial state with
%                                     maximum reach probability to identify the
%                                     vertices of the underapproximative
%                                     polytope for the stochastic reach set
%           2. init_safe_set_affine - [MUST HAVE] Affine constraints (if any) on
%                                     the initial state. Must include a
%                                     translate of the affine hull of the
%                                     set_of_dir_vecs | On intersection with the
%                                     safe set, it should result in a 2-D set
%           3. verbose              - Verbosity of the implementation {0,1}
%                                       0 - No output 
%                                       1 - Outputs the direction vector being
%                                           analyzed and the method used to
%                                           obtain xmax under study (maximizing
%                                           the reach probability or the
%                                           Chebyshev centering)
%           4. pwa_accuracy         - Accuracy of the piecewise affine
%                                     overapproximation of the inverse of the
%                                     standard normal cumulative density
%                                     function
% 
%       'genzps-open' : Genz's algorithm + Patternsearch
%           1. set_of_dir_vecs      - [MUST HAVE] Set of direction vectors shot
%                                     outwards from the initial state with
%                                     maximum reach probability to identify the
%                                     vertices of the underapproximative
%                                     polytope for the stochastic reach set
%           2. init_safe_set_affine - [MUST HAVE] Affine constraints (if any) on
%                                     the initial state. Must include a
%                                     translate of the affine hull of the
%                                     set_of_dir_vecs | On intersection with the
%                                     safe set, it should result in a 2-D set
%           3. verbose              - Verbosity of the implementation {0,1}
%                                       0 - No output 
%                                       1 - Outputs the direction vector being
%                                           analyzed, the summarized progress
%                                           made by the bisection used for line
%                                           search
%           4. tol_bisect           - Tolerance for the bisection to terminate
%                                     the line search along a direction vector
%                                     for the vertex of the underapproximative
%                                     polytope [Default 1e-2]
%           5. desired_accuracy     - Accuracy expected for the integral of the
%                                     Gaussian random vector X over the safety
%                                     tube => Accuracy of the result [Default
%                                     1e-3]
%           6. PSoptions            - MATLAB struct from psoptimset(), options
%                                     for MATLAB's patternsearch 
%                                     [Default psoptimset('Display', 'off')]
% 
%       'lag-over'/'lag-under' : Lagrangian-based over- and underapproximation
%            1. bound_set_method        - Method for obtaining the bounded set
%                                         for over or underapproximation. The
%                                         available methods are: 'polytope',
%                                         'ellipsoid', and 'load', some of which
%                                         requires additional arguments.
%               a. polytope             - Scale a user-provided polytope to
%                                         satisfy the given probability
%                                         constraint
%                            - template_polytope 
%                                       : [MUST HAVE] Template polytope
%                                         which is scaled to get the
%                                         bounded set
%                            - n_particles
%                                       : Number of particles to use for
%                                         RandomVector/getProbPolyhedron
%               b. ellipsoid            - Construct an ellipsoid to satify the
%                                         given probability constraint
%               c. load                 - Load a predefined polyhedron bounding
%                                         set; primarily used for comparison and
%                                         repeatability testing.
%                            - load_str : [MUST HAVE] Path to the file to load.
%                                         All other inputs are IRRELEVANT for
%                                         this option.  Mat files to be loaded
%                                         must have only contain the bounded set
%                                         (Polyhedron object).
%            2. verbose                 - Verbosity of the implementation
%                                         {0,1,2,3}
%                                         0 - No output 
%                                         1 - Provides feedback on the progress
%                                             of the direction vector
%                                             computation
%                                         2 - Provides timing information about
%                                             each step
%                                         3 - Provides plots of the intermediate
%                                             recursions
%            3. compute_style           - Computation style for the
%                                         set-operation methods
%               a. 'VHmethod'           - [DEFAULT] Use MPT3's Polyhedron 
%                                         manipulations to implement the 
%                                         set operations-based recursion. This
%                                         approach will fast in low dimensions,
%                                         it does not scale well with dimension
%                                         due to the vertex-facet enumeration.
%               b. 'support'            - Use support functions and convex
%                                         optimization to perform the
%                                         computations
%                            - system   : [MUST HAVE] LtvSystem/LtiSystem object
%                                         that is being analyzed
%                            - equi_dir_vecs
%                                       : [Auto-generated] Directions that
%                                         are used for vertex representation 
%                                         creation
%                            - n_underapprox_vertices 
%                                       : [Required only for lag-under] Number
%                                         of vertices to use when
%                                         underapproximating one-step backward
%                                         reach set computation for lag-under
%                                         approach for scalability
%
% Outputs:
% --------
%   options     - Collection of user-specified options
%                 (Matlab struct created using SReachSetOptions)
%
% See also SReachSet.
%
% Notes:
% * Requires init_safe_set_affine and set_of_dir_vecs for the methods 
%   'chance-open' and 'genzps-open'.
% * Requires load_str for the method lag-under with bound_set_method 'load'
% * Requires template_polytope for the method lag-under with
%   bound_set_method 'polytope'
%       - The template polytope must contain the origin
% * compute_style governs the computation style used to implement the Lagrangian
%   technique for over and under approximation. 
%       - By default, the compute_style is 'vhmethod'. This approach will fast
%         in low dimensions, it does not scale well with dimension due to the
%         vertex-facet enumeration. It implements the recursion using the native
%         operations available via MPT3.
%       - Alternatively, the compute_style 'support' uses support function to
%         implement the Lagrangian approximations. 
%           * Since polytopes are obtained by sampling the support function, we
%             require vectors that are equally spaced apart in the required
%             Euclidean space.
%               - In SReachSetOptions, these vectors are stored in the field
%               `equi_dir_vecs`.
%               - They are auto-generated by the options difference-of-convex
%                 programming.
%           * Requires system for the computation of these vectors.
%           * For lag-over, this computation style is recursion-free. 
%               - Utilizes the support function available as a simple LP
%               - Requires equally spaced vectors in R^(sys.state_dim), obtained
%                 via difference-of-convex programming.
%               - Computes an affine transformation of these vectors for
%                 meaningful overapproximation. A maximum volume ellipsoid
%                 approximately inscribed within the polytope is computed using
%                 scenario-based robust convex programming.
%               - Returns a tight overapproximation of the overapproximation
%                 polytope by sampling this support function.
%           * For lag-under, this computation style is `vhmethod guided by
%             support functions` to reduce some of the computational overhead
%             associated with `vhmethod`.
%               - The polytopes are always expressed in the half-space form with
%                 the Minowksi sum computed by projecting the higher dimensional
%                 polytope constructed in state_space x input_space into the
%                 state_space. 
%               - The half-space form of the high dimension polytope is sampled
%                 to obtain the vertex form of an underapproximate polytope,
%                 which after projection, is converted back into its half-space
%                 form.
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
    valid_bound_method = {'load','polytope','ellipsoid'};
    valid_compute_style = {'VHmethod','support'};
    % Input parsing
    inpar = inputParser();
    inpar.addRequired('prob_str', @(x) any(validatestring(x,valid_prob)));
    inpar.addRequired('method_str', @(x) any(validatestring(x,valid_method)));

    if nargin >= 2 && contains(method_str,'lag-')
        inpar.addParameter('bound_set_method',[], @(x) any(validatestring(x, ...
            valid_bound_method)));
        %% Optional arguments that are made REQUIRED based on bound_set_method
        % Get load string for load option
        inpar.addParameter('load_str', ' ', @(x) validateattributes(x, ...
            {'char'}, {'nonempty'}));
        % Equi-directional vector generation
        inpar.addParameter('n_underapprox_vertices', 100,...
            @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer',...
            'positive'}));
        % Verbosity
        inpar.addParameter('verbose',0, @(x) validateattributes(x, ...
            {'numeric'}, {'scalar','>=',0,'<=',3}));
        % System object for computation of equi-directional vector
        % generation
        inpar.addParameter('system', [],...
            @(x) validateattributes(x,{'LtvSystem','LtiSystem'}, {'nonempty'}));
        % System object for computation of equi-directional vector
        % generation
        inpar.addParameter('equi_dir_vecs', [],...
            @(x) validateattributes(x,{'numeric'}, {'nonempty'}));
        % Template polytope
        inpar.addParameter('template_polytope',[], @(x) validateattributes(x,...
            {'Polyhedron'}, {'nonempty'}));
        % n_particles
        inpar.addParameter('n_particles',1e6, @(x) validateattributes(x, ...
            {'numeric'}, {'scalar','integer','nonempty'}));
        % compute_style
        inpar.addParameter('compute_style', 'vhmethod',...
            @(x) any(validatestring(x,valid_compute_style)));
    elseif nargin >= 2
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
    else
        throwAsCaller(SrtInvalidArgsError('Too few arguments given'));
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
        if strcmpi(options.compute_style,'support') 
            if any(contains(inpar.UsingDefaults,'system'))
                throw(SrtInvalidArgsError(['system ',...
                    '(LtvSystem/LtiSystem object) is a required ',...
                    'input for SReachSet when using ',...
                    '''lag-over''/''lag-under''.']));
            end
            % get underapproximated level set (robust effective target)
            if options.verbose >= 2
                timerVal = tic;
            end
            options.equi_dir_vecs = spreadPointsOnUnitSphere(...,
                options.system.state_dim + options.system.input_dim,...
                options.n_underapprox_vertices, options.verbose);                
            if options.verbose >= 2
                fprintf('Time to spread the vectors: %1.3f s\n\n',...
                    toc(timerVal));
            end                    
        end
        switch(lower(options.bound_set_method))
            case 'load'
                if any(contains(inpar.UsingDefaults,'load_str'))
                    throw(SrtInvalidArgsError(['Expected a path to ', ...
                        'matfile with bounded set (only) as input for ', ...
                        'bound_set_method: load']));
                end
            case 'polytope'
                if any(contains(inpar.UsingDefaults,'template_polytope'))
                    throw(SrtInvalidArgsError(['Expected template_polytope,',...
                        ' a Polyhedron object as an input for ', ...
                        'bound_set_method: polytope']));
                end
                if ~options.template_polytope.contains(...
                        zeros(options.template_polytope.Dim,1))
                    throw(SrtInvalidArgsError(['Expected template_polytope,',...
                        ' a Polyhedron object that contains 0 as input for ', ...
                        'bound_set_method: polytope']));
                end
            case 'ellipsoid'        
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
