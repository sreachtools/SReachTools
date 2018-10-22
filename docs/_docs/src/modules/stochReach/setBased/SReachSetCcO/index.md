---
layout: docs
title: SReachSetCcO.m
---

```
  Obtain an open-loop controller-based underaproximative stochastic reach-avoid
  set using Fourier transform, convex optimization, and patternsearch
  =============================================================================
 
  getUnderapproxStochReachAvoidSet computes the open-loop controller-based
  underapproximative stochastic reach-avoid set  to the terminal hitting-time
  stochastic reach-avoid problem discussed in
 
  A. Vinod, and M. Oishi, "Scalable Underapproximative Verification of
  Stochastic LTI Systems Using Convexity and Compactness," in Proceedings of
  Hybrid Systems: Computation and Control (HSCC), 2018. 
 
  Specifically, we use chance-constrained approach to speed up the computation.
 
  USAGE: TODO
 
  =============================================================================
 
  [underapprox_stoch_reach_avoid_polytope, ...
   opt_input_vector_at_vertices, ...
   varargout] = ...
   getUnderapproxStochReachAvoidSet(...
                                 sys, ...
                                 safety_tube, ...
                                 init_safe_set_affine_const, ...
                                 prob_thresh, ...
                                 set_of_dir_vecs, ...
                                 varargin)
  
  Inputs:
  -------
    sys                  - LtiSystem object describing the system to be verified
    safety_tube          - Target tube to stay within [Tube object]
    init_safe_set_affine_const        
                         - Affine constraints (if any) on the initial state
                           Must include a translate of the affine hull of the
                           set_of_dir_vecs                          
    prob_thresh 
                         - Probability thresh (\theta) that defines the
                           stochastic reach-avoid set 
                           {x_0: V_0^\ast( x_0) \geq \theta}
    set_of_dir_vecs
                         - Number of unique directions defining the polytope
                           vertices. Its span is the affine hull whose slice of
                           the stochastic reach-avoid set is of interest.
    method               - TODO
    options.desired_accuracy     - (Optional for 'genzps') Accuracy expected for the
                           integral of the Gaussian random vector X over the
                           concatenated_safety_tube [Default 5e-3]
    PSoptions            - (Optional for 'genzps') Options for patternsearch 
                           [Default psoptimset('Display', 'off')]
 
  Outputs:
  --------
    underapprox_stoch_reach_avoid_polytope
                         - Underapproximative polytope of dimension
                           sys.state_dim which underapproximates the
                           terminal-hitting stochastic reach avoid set
    opt_input_vector_at_vertices 
                         - Optimal open-loop policy ((sys.input_dim) *
                           time_horizon)-dim.  vector U = [u_0; u_1; ...; u_N]
                           (column vector) for each vertex of the polytope
    xmax                 - (Optional) Initial state that has the maximum
                           stochastic reach-avoid prob using an open-loop
                           controller
    opt_input_vector_for_xmax
                         - (Optional) Optimal open-loop policy
                           ((sys.input_dim) * time_horizon)-dimensional
                           vector U = [u_0; u_1; ...; u_N] (column vector) for
                           xmax
    max_underapprox_reach_avoid_prob
                         - (Optional) Maximum attainable stochastic reach-avoid
                           prob using an open-loop controller; Maximum
                           terminal-hitting time reach-avoid prob at xmax
    opt_theta_i          - (Optional) Vector comprising of scaling factors along
                           each direction of interest
    opt_reachAvoid_i     - (Optional) Maximum terminal-hitting time reach-avoid
                           prob at the vertices of the polytope
    vertex_underapprox_polytope
                         - (Optional) Vertices of the polytope: xmax + opt_theta_i *
                           set_of_dir_vecs
    R                    - (Optional for ccc) Chebyshev radius associated with
                           xmax
 
  See also examples/FtCVXUnderapproxVerifyCWH.mlx*.
 
  Notes:
  ------
  * NOT ACTIVELY TESTED: Builds on other tested functions.
  * MATLAB DEPENDENCY: Uses MATLAB's Global Optimization Toolbox; Statistics and
                       Machine Learning Toolbox.
                       Needs patternsearch for gradient-free optimization
                       Needs normpdf, normcdf, norminv for Genz's algorithm
  * EXTERNAL DEPENDENCY: Uses MPT3 and CVX
                       Needs MPT3 for defining a controlled system and the
                       definition of the safe, the target (polytopic) sets, and
                       the affine hull of interest
                       Needs CVX to setup convex optimization problems that
                       1) initializes the patternsearch-based optimization, and
                       2) computes the upper bound for the bisection
  * Specify both options.desired_accuracy and PSoptions or neither to use the defaults 
  * max_underapprox_reach_avoid_prob is the highest thresh
    that may be given while obtaining a non-trivial underapproximation
  * See @LtiSystem/getConcatMats for more information about the
      notation used.
  
  =============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
 
```
