---
layout: docs
title: SReachSetOptions.m
---

```
  Create user-specifiable options for use with SReachSet()
  =============================================================================
 
  SReachSetOptions creates a MATLAB struct that contains user-specifiable
  options that may be used with SReachSet
 
  =============================================================================
 
    options = SReachSetOptions(prob_str, method_str, varargin)
  
  Inputs:
  -------
    prob_str    - String specifying the problem of interest. For each case, we
                  compute the optimal value function that maps initial states
                  to different maximal reach probabilities
                      1. 'term' : Stay within the safety_tube
    method_str  - Solution technique to be used; available techniques:
                      'chance-opn', 'genzps-open', 'lag-under', 'lag-over'
    varargin    - Additional required options for each technique, specified as
                  Name-Value pairs. The additional required options are 
                  specified below.
  
        'chance-open' : Convex chance-constrained approach for an open-loop 
                        controller synthesis
            1. set_of_dir_vecs      - [MUST HAVE] Set of direction vectors shot
                                      outwards from the initial state with
                                      maximum reach probability to identify the
                                      vertices of the underapproximative
                                      polytope for the stochastic reach set
            2. init_safe_set_affine - [MUST HAVE] Affine constraints (if any) on
                                      the initial state. Must include a
                                      translate of the affine hull of the
                                      set_of_dir_vecs | On intersection with the
                                      safe set, it should result in a 2-D set
            3. verbose              - Verbosity of the implementation {0,1}
                                        0 - No output 
                                        1 - Outputs the direction vector being
                                            analyzed and the method used to
                                            obtain xmax under study (maximizing
                                            the reach probability or the
                                            Chebyshev centering)
            4. pwa_accuracy         - Accuracy of the piecewise affine
                                      overapproximation of the inverse of the
                                      standard normal cumulative density
                                      function
  
        'genzps-open' : Genz's algorithm + Patternsearch
            1. set_of_dir_vecs      - [MUST HAVE] Set of direction vectors shot
                                      outwards from the initial state with
                                      maximum reach probability to identify the
                                      vertices of the underapproximative
                                      polytope for the stochastic reach set
            2. init_safe_set_affine - [MUST HAVE] Affine constraints (if any) on
                                      the initial state. Must include a
                                      translate of the affine hull of the
                                      set_of_dir_vecs | On intersection with the
                                      safe set, it should result in a 2-D set
            3. verbose              - Verbosity of the implementation {0,1}
                                        0 - No output 
                                        1 - Outputs the direction vector being
                                            analyzed, the summarized progress
                                            made by the bisection used for line
                                            search
            4. tol_bisect           - Tolerance for the bisection to terminate
                                      the line search along a direction vector
                                      for the vertex of the underapproximative
                                      polytope [Default 1e-2]
            5. desired_accuracy     - Accuracy expected for the integral of the
                                      Gaussian random vector X over the safety
                                      tube => Accuracy of the result [Default
                                      1e-3]
            6. PSoptions            - MATLAB struct from psoptimset(), options
                                      for MATLAB's patternsearch 
                                      [Default psoptimset('Display', 'off')]
  
        'lag-over'/'lag-under' : Lagrangian-based over- and underapproximation
             1. bound_set_method        - Method for obtaining the bounded set
                                          for over or underapproximation. The
                                          available methods are: 'polytope',
                                          'ellipsoid', and 'load', some of which
                                          requires additional arguments.
                a. polytope             - Scale a user-provided polytope to
                                          satisfy the given probability
                                          constraint
                             - template_polytope 
                                        : [MUST HAVE] Template polytope
                                          which is scaled to get the
                                          bounded set
                             - n_particles
                                        : Number of particles to use for
                                          RandomVector/getProbPolyhedron
                b. ellipsoid            - Construct an ellipsoid to satify the
                                          given probability constraint
                             - system   : [MUST HAVE] LtvSystem/LtiSystem object
                                          that is being analyzed
                             - n_underapprox_vertices 
                                        : Number of vertices to use when
                                          underapproximating one-step backward
                                          reach set computation for lag-under
                                          approach for scalability
                             - equi_dir_vecs
                                        : [Auto-generated] Directions that
                                          are used for vertex representation 
                                          creation
                c. load                 - Load a predefined polyhedron bounding
                                          set; primarily used for comparison and
                                          repeatability testing.
                             - load_str : [MUST HAVE] Path to the file to load.
                                          All other inputs are IRRELEVANT for
                                          this option.  Mat files to be loaded
                                          must have only contain the bounded set
                                          (Polyhedron object).
             2. verbose                 - Verbosity of the implementation
                                          {0,1,2,3}
                                          0 - No output 
                                          1 - Provides feedback on the progress
                                              of the direction vector
                                              computation
                                          2 - Provides timing information about
                                              each step
                                          3 - Provides plots of the intermediate
                                              recursions
 
  Outputs:
  --------
    options     - Collection of user-specified options
                  (Matlab struct created using SReachSetOptions)
 
  See also SReachSet.
 
  Notes:
  * Requires init_safe_set_affine and set_of_dir_vecs for the methods 
    'chance-open' and 'genzps-open'.
  * Requires load_str for the method lag-under with bound_set_method 'load'
  * Requires template_polytope for the method lag-under with
    bound_set_method 'polytope'
        - The template polytope must contain the origin
  * Requires system for the method lag-under with bound_set_method 'ellipsoid'
  * equi_dir_vecs are auto-generated only if the bound_set_method is 'ellipsoid'
  
  ============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```
