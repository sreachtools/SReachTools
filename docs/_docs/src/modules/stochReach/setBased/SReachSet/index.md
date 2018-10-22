---
layout: docs
title: SReachSet.m
---

```
  Compute the stochastic reach set corresponding to the stochastic reachability 
  problem of a target tube using a host of techniques
  =============================================================================
 
  SReachPoint computes an approximation to the stochastic reach set of the 
  stochastic reachability problem of a target tube. Specifically, it
  computes an approximation to the set {x_0 \in X: V_0(x_0) >= \alpha}
  where
 
      V_0(x_0) = maximize Prob( \cap_{i=1}^N x_t lies in Safe_t)
                 subject to
                    dynamics and bounds on control
 
  In other words, V_0 is the optimal value function corresponding to the 
  terminal hitting-time stochastic reachability problem (stochastic
  reachability of a target tube).
 
  We use the theory discussed in,
  
  1. A. P. Vinod and M. Oishi, HSCC 2018 TODO
  2. J. Gleason, A. P. Vinod, and M. Oishi, CDC 2017 TODO
  
  to compute these sets
 
 
  This function is a compilation of various techniques proposed in the
  literature:
 
  1. Convex chance-constrained-based approach (chance-open):
 
     High-level desc.   : Use Boole's inequality, Gaussian random vector, and
                          piecewise linear approximation of the inverse of the
                          standard normal cumulative density function to create
                          a linear program-based approximation to the original
                          optimization
     Approximation      : Guaranteed underapproximation
     Controller type    : Open-loop controller that satisfies the hard
                          input bounds 
     SReachTool function: SReachSetCcO
     Dependency (EXT)   : CVX
     Dependency (MATLAB): Symbolic toolbox
     Paper              : a. Lesser, Oishi, Erwin TODO.
                          b. A. Vinod and M. Oishi, HSCC 2018 TODO
 
  2. Fourier transform + Patternsearch (genzps-open):
 
     High-level desc.   : Maximize the multivariate Gaussian integral over a
                          polytope, evaluated using Genz's algorithm, and
                          optimize the nonlinear (log-concave) problem using
                          MATLAB's patternsearch
     Approximation      : Approximate upto a user-specified tolerance
     Controller type    : Open-loop controller that satisfies the hard input
                          bounds
     Dependency (MATLAB): Global Optimization toolbox (for patternsearch)
     SReachTool function: SReachSetGpO
     Paper              : b. A. Vinod and M. Oishi, HSCC 2018 TODO
 
  3. Lagrangian underapproximation
 
     High-level desc.   : Use computational geometry tools to compute a
                          robust effective target set
     Approximation      : Guaranteed underapproximation
     Controller type    : Closed-loop controller that satisfies the hard input
                          bounds
     SReachTool function: SReachSetLag
     Paper              : J. Gleason, A. P. Vinod, and M. Oishi, CDC 2017 TODO
 
  4. Lagrangian overapproximation
 
     High-level desc.   : Use computational geometry tools to compute a
                          augmented effective target set
     Approximation      : Guaranteed overapproximation
     Controller type    : Closed-loop controller that satisfies the hard input
                          bounds
     SReachTool function: SReachSetLag
     Paper              : J. Gleason, A. P. Vinod, and M. Oishi, Automatica TODO
 
 
  USAGE: TODO
 
  =============================================================================
 
  [stoch_reach_prob, opt_controller, varargout] = SReachSet(prob_str, ...
     method_str, sys, initial_state, safety_tube, options)
  
  Inputs:
  -------
    prob_str    - String specifying the problem of interest. For each case, we
                  compute the optimal value function that maps initial states
                  to different maximal reach probabilities
                      1. 'term' : Stay within the safety_tube
    method_str  - Solution technique to be used.
                      'chance-open' -- Convex chance-constrained approach for an
                                       open-loop controller synthesis
                      'genzps-open' -- Genz's algorithm + Patternsearch
                      'lag-under'   -- Set-computation (Lagrangian)-based
                                       underapproximation
                      'lag-over'    -- Set-computation (Lagrangian)-based 
                                       overapproximation
    sys         - System description (LtvSystem/LtiSystem object)
    prob_thresh - Probability threshold at which the set is to be constructed
    safety_tube - Collection of (potentially time-varying) safe sets that
                  define the safe states (Tube object)
    options     - Collection of user-specified options for each of the solution
                  (Matlab struct created using SReachSetOptions)
 
  Outputs:
  --------
    stoch_reach_set 
                - Approximation (underapproximation, in some cases) of the
                  stochastic reach set
    other_info_polytope
                - [Available for 'chance-open'/'genzps-open'] A MATLAB struct
                  containing additional info:
                  1. opt_input_vec_vertex 
                            Optimal open-loop policy ((sys.input_dim) *
                            time_horizon)-dim. vector U = [u_0; u_1; ...; u_N]
                            (column vector) for each vertex of the polytope
                  2. xmax - Initial state that has the maximum stochastic reach
                            probability using an open-loop controller (via the
                            method in use)
                  3. opt_input_vec_xmax
                          - Optimal open-loop policy ((sys.input_dim) *
                            time_horizon)-dimensional vector U = [u_0; u_1;
                            ...; u_N] (column vector) for xmax (via the method
                            in use)
                  4. opt_underapprox_reach_prob_xmax
                          - Maximum attainable stochastic reach probability
                            using an open-loop controller (via the method in
                            use)
                  5. opt_theta_i 
                          - Vector comprising of scaling factors along each
                            direction of interest 
                  6. opt_reach_prob_i     
                          - Maximum terminal-hitting time reach probability at
                            the vertices of the polytope
                  7. vertex_i    
                          - Vertices of the polytope: xmax + opt_theta_i *
                            set_of_dir_vecs
                  8. R    - Chebyshev radius associated with xmax
    other_info_lag
                - [Available for 'lag-X'] TODO
 
  Notes:
  * See @LtiSystem/getConcatMats for more information about the notation used.
  
  =============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
 
```