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
                      'chance-open', 'genzps-open', 'lag-under', 'lag-over'
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
                                      A (sys.state_dim x n_dir)-dimensional
                                      matrix
            2. init_safe_set_affine - [MUST HAVE] Affine constraints (if any) on
                                      the initial state. Must include a
                                      translate of the affine hull of the
                                      set_of_dir_vecs | On intersection with the
                                      safe set at t=0, it should result in a 2-D 
                                      set
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
                                      5e-2] | This value can't be smaller
                                      than 1e-2
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
                             - desired_accuracy
                                        : Accuracy for 
                                          RandomVector/getProbPolyhedron | This 
                                          value can't be smaller than 1e-2
                                          [Default 1e-2]
                b. ellipsoid            - Construct an ellipsoid to satify the
                                          given probability constraint
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
             3. compute_style           - Computation style for the
                                          set-operation methods
                a. 'vfmethod'          - [DEFAULT] Use MPT3's Polyhedron 
                                          manipulations to implement the 
                                          set operations-based recursion. This
                                          approach will fast in low dimensions,
                                          it does not scale well with dimension
                                          due to the vertex-facet enumeration.
                b. 'support'            - Use support functions and convex
                                          optimization to perform the
                                          computations
                             - system   : [MUST HAVE] LtvSystem/LtiSystem object
                                          that is being analyzed
                             - n_vertices
                                        : Number of vertices to use 
                                             [For lag-under] underapproximating
                                                    one-step backward reach set
                                             [For lag-under] overapproximating
                                                    the stochastic reach set
                                                    overapproximation 
                             - equi_dir_vecs
                                        : [Auto-generated] Directions that
                                          are used for overapproximation
             4. vf_enum_method          - Enumeration style to use for
                                          vertex-facet enumeration. See notes
                a. 'cdd'                - [DEFAULT] Use MPT3's native polyhedral
                                          operations, which in turn use CDDMEX
                                          for vertex-facet enumeration.
                b. 'lrs'                - Use GeoCalcLib (a MATLAB bridge to
                                          McGill's LRS C code base) for
                                          vertex-facet enumeration
 
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
  * compute_style governs the computation style used to implement the Lagrangian
    technique for over and under approximation. 
        - By default, the compute_style is 'vfmethod'. This approach will fast
          in low dimensions, it does not scale well with dimension due to the
          vertex-facet enumeration. It implements the recursion using the native
          operations available via MPT3.
        - Alternatively, the compute_style 'support' uses support function to
          implement the Lagrangian approximations. 
            * Since polytopes are obtained by sampling the support function, we
              require vectors that are equally spaced apart in the required
              Euclidean space.
                - In SReachSetOptions, these vectors are stored in the field
                `equi_dir_vecs`.
                - They are auto-generated by the options difference-of-convex
                  programming.
            * Requires system for the computation of these vectors.
            * For lag-over, this computation style is recursion-free. 
                - Utilizes the support function available as a simple LP
                - Requires equally spaced vectors in R^(sys.state_dim), obtained
                  via difference-of-convex programming.
                - Computes an affine transformation of these vectors for
                  meaningful overapproximation. A maximum volume ellipsoid
                  approximately inscribed within the polytope is computed using
                  scenario-based robust convex programming.
                - Returns a tight overapproximation of the overapproximation
                  polytope by sampling this support function.
            * For lag-under, this computation style is `vfmethod guided by
              support functions` to reduce some of the computational overhead
              associated with `vfmethod`.
                - The polytopes are always expressed in the half-space form with
                  the Minowksi sum computed by projecting the higher dimensional
                  polytope constructed in state_space x input_space into the
                  state_space. 
                - The half-space form of the high dimension polytope is sampled
                  to obtain the vertex form of an underapproximate polytope,
                  which after projection, is converted back into its half-space
                  form.
  * While specifying n_vertices, use the formula 
    2^{n_dim} points_per_quad + 2 * n_dim 
    to obtain a spread of points where each quadrant has `points_per_quad`
    and the standard axis are also included.
        - For lag-under, n_dim is sys.state_dim + sys.input_dim
        - For lag-over, n_dim is sys.state_dim
  * Vertex-facet enumeration is a computationally hard problem. SReachTools
    currently support two popular techniques for addressing the same:
    1. CDDMEX - MPT's preferred approach for vertex-facet enumeration
                * Requires no additional installation steps
                * Known to provide incorrect results or fail completely in some
                  cases
                * See following websites for more information: 
                    https://www.inf.ethz.ch/personal/fukudak/cdd_home/index.html
                    http://www.swmath.org/software/5097
    2. LRS    - Avis's LRS with MATLAB interface provided Rainer's GeoCalcLib
                * Requires few additional installation steps
                * Worked more reliably than CDDMEX
                * See following websites for more information: 
                    http://cgm.cs.mcgill.ca/~avis/C/lrs.html
                    http://worc4021.github.io/GeoCalcLib/
  * While 'Gaussian' disturbance can have options.bound_set_method be 'polytope'
    or 'ellipsoid', 'UserDefined' disturbance requires options.bound_set_method
    to be 'polytope'.
  ============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```
