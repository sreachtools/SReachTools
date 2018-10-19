---
layout: docs
title: SReachPointPaO.m
---

```
  Solve the stochastic reach-avoid problem (lower bound on the probability and
  open-loop controller synthesis) using particle filter control
  =============================================================================
 
  SReachPointCcO implements a mixed-integer linear program-based approximation
  to the stochastic reachability of a target tube problem. This solution is
  based off the particle filter control formulation discussed in
 
  K. Lesser, M. Oishi, and R. Erwin, "Stochastic reachability for control of
  spacecraft relative motion," in IEEE Conference on Decision and Control (CDC),
  2013.
 
  =============================================================================
 
  [approx_stoch_reach, opt_input_vec, risk_alloc_state, varargout] =...
     SReachPointCcO(sys, initial_state, safety_tube, options)
 
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
                    tube problem computed using particle control
    opt_input_vec
                - Open-loop controller: column vector of dimension
                  (sys.input_dim*N) x 1
 
  See also SReachPoint.
 
  Notes:
  * Requires Gurobi as the backend solver for optimizing the resulting
        mixed-integer linear program
  * See @LtiSystem/getConcatMats for more information about the notation used.
  
  ============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```
