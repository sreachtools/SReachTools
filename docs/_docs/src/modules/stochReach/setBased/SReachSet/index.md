---
layout: docs
title: SReachSet.m
---

```
  Compute the stochastic reach set corresponding to the stochastic reachability 
  problem of a target tube using a host of techniques
  =============================================================================
 
  SReachSet computes an approximation to the stochastic reach set of the 
  stochastic reachability problem of a target tube. Specifically, it
  computes an (over- and under-)approximation to the set 
   
     L_0(\alpha) = {x_0 \in X: V_0(x_0) >= \alpha} 
 
  where
 
      V_0(x_0) = maximize Prob( \cap_{i=1}^N x_t lies in Safe_t)
                 subject to
                    dynamics and bounds on control
 
  In other words, V_0 is the maximal reach probability associated with the
  problem of stochastic reachability of a target tube, and L_0(\alpha) is the
  stochastic reach set --- the set of initial states from which an admissible
  controller exists that can drive the state to stay within the safety tube
  {Safe_t}_{t=1}^{N}.
 
  This function is a compilation of various techniques proposed in the
  literature:
 
  1. Convex chance-constrained-based approach (chance-open):
 
     Set computation    : Perform line searches along a set of user-specified
                          direction vectors originating from the best performing
                          initial state; Point-based stochastic reachability via
                          chance-constraint formulation
     Approximation      : Guaranteed underapproximation
     Paper              : a. A. Vinod and M. Oishi, "Scalable underapproximative
                             verification of stochastic LTI systems using
                             convexity and compactness," In Proc. Hybrid Syst.:
                             Comput. & Ctrl., pages 1--10, 2018. HSCC 2018
                          b. A. Vinod and M. Oishi, "Stochastic reachability of
                             a target tube: Theory and computation," IEEE
                             Transactions in Automatic Control, 2018 (submitted)
                             https://arxiv.org/pdf/1810.05217.pdf.
     See also SReachPointCcO.
 
  2. Fourier transform + Patternsearch (genzps-open):
 
     Set computation    : Perform line searches along a set of user-specified
                          direction vectors originating from the best performing
                          initial state; Point-based stochastic reachability via
                          optimization of the multivariate Gaussian integral
                          over a polytope (Genz's algorithm) using MATLAB's
                          patternsearch
     Approximation      : Approximate upto a user-specified tolerance
     Controller type    : Open-loop controller that satisfies the hard input
                          bounds
     Dependency (MATLAB): Global Optimization toolbox (for patternsearch)
     SReachTool function: SReachSetGpO
     Paper              : a. A. Vinod and M. Oishi, "Scalable underapproximative
                             verification of stochastic LTI systems using
                             convexity and compactness," In Proc. Hybrid Syst.:
                             Comput. & Ctrl., pages 1--10, 2018. HSCC 2018
                          b. A. Vinod and M. Oishi, "Scalable Underapproximation
                             for Stochastic Reach-Avoid Problem for
                             High-Dimensional LTI Systems using Fourier
                             Transforms," in IEEE Control Systems Letters
                             (L-CSS), 2017.
 
  3. Lagrangian underapproximation
 
     High-level desc.   : Use computational geometry tools to compute an
                          underapproximation of the stochastic reach set
     Approximation      : Guaranteed underapproximation
     Controller type    : Closed-loop controller that satisfies the hard input
                          bounds
     SReachTool function: SReachSetLag
     Paper              : a. J. Gleason, A. Vinod, and M. Oishi, "Lagrangian
                             Approximations for Stochastic Reachability of a
                             Target Tube," 2018.
                             https://arxiv.org/abs/1810.07118 TODO
                          b. J. Gleason, A. Vinod, and M. Oishi,
                             "Underapproximation of Reach-Avoid Sets for
                             Discrete-Time Stochastic Systems via Lagrangian
                             Methods," In Proceedings of the IEEE Conference on
                             Decision and Control, 2017
 
  4. Lagrangian overapproximation
 
     High-level desc.   : Use computational geometry tools to compute an
                          overapproximation of the stochastic reach set
     Approximation      : Guaranteed overapproximation
     Controller type    : Closed-loop controller that satisfies the hard input
                          bounds
     SReachTool function: SReachSetLag
     Paper              : a. J. Gleason, A. Vinod, and M. Oishi, "Lagrangian
                             Approximations for Stochastic Reachability of a
                             Target Tube," 2018.
                             https://arxiv.org/abs/1810.07118
                          b. J. Gleason, A. Vinod, and M. Oishi,
                             "Underapproximation of Reach-Avoid Sets for
                             Discrete-Time Stochastic Systems via Lagrangian
                             Methods," In Proceedings of the IEEE Conference on
                             Decision and Control, 2017
 
  See also examples/cwhSReachSetDemo.m and examples/dubinsSReachSetDemo.m.
 
  =============================================================================
 
  [stoch_reach_set, [extra_info]] = SReachSet(prob_str, method_str, sys, ...
        prob_thresh, safety_tube, options)
  
  Inputs:
  -------
    prob_str    - String specifying the problem of interest. For each case, we
                  compute the optimal value function that maps initial states
                  to different maximal reach probabilities
                      1. 'term' : Stay within the safety_tube
    method_str  - Solution technique to be used.
                      'chance-open' -- Underapproximative construction of the
                                       stochastic reach set using convex
                                       chance-constrained approach
                      'genzps-open' -- Underapproximative construction of the
                                       stochastic reach set using Genz's
                                       algorithm + Patternsearch |
                                       Underapproximation holds true upto a
                                       user-specified error for the multivariate
                                       integral
                      'lag-under'   -- Underapproximative construction of the
                                       stochastic reach set using set-theoretic
                                       (Lagrangian) approach
                      'lag-over'    -- Overapproximative construction of the
                                       stochastic reach set using set-theoretic
                                       (Lagrangian) approach
    sys         - System description (LtvSystem/LtiSystem object)
    prob_thresh - Probability threshold at which the set is to be constructed
    safety_tube - Collection of (potentially time-varying) safe sets that
                  define the safe states (Tube object)
    options     - Collection of user-specified options for each of the solution
                  (Matlab struct created using SReachSetOptions)
 
  Outputs:
  --------
    stoch_reach_set 
                - Approximation (over- or under-approximation) of the
                  stochastic reach set
    extra_info  - A MATLAB struct containing additional info, like optimal
                  open-loop input vector from the vertices and the initial state
                  with maximum reach probability in case of
                  'chance-open'/'genzps-open', and the effective_target_tube and
                  the bounded disturbance set in case of 'lag-over/lag-under'.
                  See the docstring of SReachSetXXX for more details.
 
  Notes:
  * 'set_of_dirs' and 'init_safe_set_affine' needs to be provided to the options
    if 'chance-open' or 'genzps-open' is to be used. See SReachSetOptions() for
    more details
  * See @LtiSystem/getConcatMats for more information about the notation used.
  * For lagrangian underapproximation approach, see getSReachLagUnderapprox.
      - From computational geometry, intersections and Minkowski differences are
        best performed in facet representation and Minkowski sums are best
        performed in vertex representation. However, since in this computation,
        all three operations are required, scalability of the algorithm is
        severly hampered, despite theoretical elegance.
      - Since box and random approaches in SReachSetOptions produce Polyhedron
        objects for disturbance sets, we rely on MPT for all the set operations.
        This means we do have scalability issues mentioned above.
      - For ellipsoid approach in SReachSetOptions, we seek a purely facet-based
        operation and utilize the ray-shooting algorithm to compute a
        facet-based underapproximation of the Minkowski sum step (via
        vertex-based underapproximation, followed by projection, followed by
        convex hull operation)
      - While 'Gaussian' disturbance can have options.bound_set_method be
        'polytope' or 'ellipsoid', 'UserDefined' disturbance requires
        options.bound_set_method to be 'polytope'.
 
  =============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
 
```
