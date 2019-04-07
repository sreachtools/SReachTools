---
layout: docs
title: SReachPointCcO.m
---

```
  Solve the problem of stochastic reachability of a target tube (a lower bound
  on the maximal reach probability and an open-loop controller synthesis) using
  convex chance-constrained optimization
  =============================================================================
 
  SReachPointCcO implements a chance-constrained convex underapproximation to
  the stochastic reachability of a target tube problem. The original problem was
  formulated (for the simpler problem of terminal hitting-time stochastic
  reach-avoid problem) in
 
  K. Lesser, M. Oishi, and R. Erwin, "Stochastic reachability for control of
  spacecraft relative motion," in IEEE Conference on Decision and Control (CDC),
  2013.
 
  This function implements a convex solver-friendly using piecewise-affine
  overapproximations of the convex constraints, as discussed in
 
  A. Vinod and M. Oishi. Affine controller synthesis for stochastic reachability
  via difference of convex programming. In Proc. Conf. Dec. & Ctrl., 2019.
  (submitted). https://hscl.unm.edu/affinecontrollersynthesis/
 
     High-level desc.   : Use Boole's inequality, Gaussian random vector, and
                          piecewise linear approximation of the inverse of the
                          standard normal cumulative density function to create
                          a linear program-based approximation to the original
                          optimization
     Approximation      : Guaranteed underapproximation
     Controller type    : Open-loop controller that satisfies the hard
                          input bounds 
     Optimality         : Optimal open-loop controller for the
                          underapproximation problem due to convexity guarantees
 
  =============================================================================
 
  [lb_stoch_reach, opt_input_vec, risk_alloc_state, varargout] = ...
     SReachPointCcO(sys, initial_state, safety_tube, options)
 
  Inputs:
  -------
    sys          - System description (LtvSystem/LtiSystem object)
    initial_state- Initial state for which the maximal reach probability must be
                   evaluated (A numeric vector of dimension sys.state_dim)
    safety_tube  - Collection of (potentially time-varying) safe sets that
                   define the safe states (Tube object)
    options      - Collection of user-specified options for 'chance-open'
                   (Matlab struct created using SReachPointOptions)
 
  Outputs:
  --------
    lb_stoch_reach 
                - Lower bound on the stochastic reachability of a target tube
                  problem computed using convex chance constraints and
                  piecewise affine approximation. While it is expected to lie in
                  [0,1], it is set to -1 in cases where the CVX optimization
                  fails (cvx_status \neq Solved).
    opt_input_vec
                - Open-loop controller: column vector of dimension
                  (sys.input_dim*N) x 1
    risk_alloc_state 
                - Risk allocation for the state constraints
    extra_info  - [Optional] Useful information to construct the
                    reachability problem | Used by 'genzps-open' to avoid 
                    unnecessary recomputation
                  Matlab struct with members --- concat_safety_tube_A,
                    concat_safety_tube_b, concat_input_space_A,
                    concat_input_space_b, H, mean_X_sans_input,
                    cov_X_sans_input;
 
  See also SReachPoint.
 
  Notes:
  * We recommend using this function through SReachPoint.
  * This function requires CVX to work.
  * See @LtiSystem/getConcatMats for more information about the notation used.
  
  ============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```
