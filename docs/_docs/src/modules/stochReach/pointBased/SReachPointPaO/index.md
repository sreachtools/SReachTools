---
layout: docs
title: SReachPointPaO.m
---

```
  Solve the problem of stochastic reachability of a target tube (a lower bound
  on the maximal reach probability and an open-loop controller synthesis) using
  particle filter control
  =============================================================================
 
  SReachPointPaO implements a mixed-integer linear program-based approximation
  to the stochastic reachability of a target tube problem. This solution is
  based off the particle filter control formulation (for the simpler terminal
  hitting-time stochastic reach-avoid problem) discussed in
 
  K. Lesser, M. Oishi, and R. Erwin, "Stochastic reachability for control of
  spacecraft relative motion," in IEEE Conference on Decision and Control (CDC),
  2013.
 
     High-level desc.   : Sample scenarios based on the additive noise and solve
                          a mixed-integer linear program to make the maximum
                          number of scenarios satisfy the reachability objective
     Approximation      : No direct approximation guarantees. Accuracy improves
                          as the number of scenarios considered increases.
     Controller type    : Open-loop controller that satisfies the hard input
                          bounds
     Optimality         : Optimal (w.r.t scenarios drawn) open-loop controller
                          for the underapproximation problem 
 
  =============================================================================
 
  [approx_stoch_reach, opt_input_vec] = SReachPointPaO(sys, initial_state,...
    safety_tube, options)
 
  Inputs:
  -------
    sys          - System description (LtvSystem/LtiSystem object)
    initial_state- Initial state for which the maximal reach probability must be
                   evaluated (A numeric vector of dimension sys.state_dim)
    safety_tube  - Collection of (potentially time-varying) safe sets that
                   define the safe states (Tube object)
    options      - Collection of user-specified options for 'particle-open'
                   (Matlab struct created using SReachPointOptions)
 
  Outputs:
  --------
    approx_stoch_reach 
                - An approximation of the stochastic reachability of a target
                  tube problem computed using particle control. While it is
                  expected to lie in [0,1], it is set to -1 in cases where the
                  CVX optimization fails (cvx_status \not\in {Solved,
                  Inaccurate/Solved}).
    opt_input_vec
                - Open-loop controller: column vector of dimension
                  (sys.input_dim*N) x 1
 
  See also SReachPoint.
 
  Notes:
  * This function requires CVX with Gurobi as the backend solver for optimizing
    the resulting mixed-integer linear program.
  * See @LtiSystem/getConcatMats for more information about the notation used.
  
  ============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```
