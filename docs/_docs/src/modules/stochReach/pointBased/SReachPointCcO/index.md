---
layout: docs
title: SReachPointCcO.m
---

```
  Solve the stochastic reach-avoid problem (lower bound on the probability and
  open-loop controller synthesis) using chance-constrained convex optimization
  =============================================================================
 
  SReachPointCcO implements a chance-constrained convex underapproximation to
  the stochastic reachability of a target tube problem. This solution is based
  off the formulation (A) discussed in
 
  K. Lesser, M. Oishi, and R. Erwin, "Stochastic reachability for control of
  spacecraft relative motion," in IEEE Conference on Decision and Control (CDC),
  2013.
 
  This function implements a convex solver-friendly piecewise-affine restriction
  of the formulation (A), as discussed in
 
  A. Vinod and M. Oishi, HSCC 2018 TODO
 
  =============================================================================
 
  [lb_stoch_reach, opt_input_vec, risk_alloc_state, varargout] =...
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
                    piecewise affine approximation
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
  * See @LtiSystem/getConcatMats for more information about the notation used.
  
  ============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```
