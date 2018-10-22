---
layout: docs
title: SReachPointCcA.m
---

```
  Solve the stochastic reach-avoid problem (lower bound on the probability and
  an affine controller synthesis) using chance-constrained convex optimization
  =============================================================================
 
  SReachPointCcA implements the chance-constrained convex underapproximation to 
  the problem of stochastic reachability of a target tube
 
  A. Vinod and M. Oishi, HSCC, 2019 (TODO)
 
  This function uses difference-of-convex algorithm (also known as the 
  convex-concave procedure) to compute a local optima for the risk
  allocation and an associated affine controller.
 
  =============================================================================
 
    [lb_stoch_reach, opt_input_vec, opt_input_gain, ...
        risk_alloc_state, risk_alloc_input] = SReachPointCcA(sys, ...
         initial_state, safety_tube, options)
  
  Inputs:
  -------
    sys          - System description (LtvSystem/LtiSystem object)
    initial_state- Initial state for which the maximal reach probability must be
                   evaluated (A numeric vector of dimension sys.state_dim)
    safety_tube  - Collection of (potentially time-varying) safe sets that
                   define the safe states (Tube object)
    options      - Collection of user-specified options for 'chance-affine'
                   (Matlab struct created using SReachPointOptions)
 
  Outputs:
  --------
    lb_stoch_reach 
                - Lower bound on the stochastic reachability of a target 
                  tube problem computed using convex chance
                       constraints and difference-of-convex techniques
    opt_input_vec, 
      opt_input_gain
                - Controller U=MW+d for a concatenated input vector 
                    U = [u_0; u_1; ...; u_{N-1}] and concatenated disturbance
                    vector W=[w_0; w_1; ...; w_{N-1}]. 
                    - opt_input_gain: Affine controller gain matrix of dimension
                        (sys.input_dim*N) x (sys.dist_dim*N)
                    - opt_input_vec: Open-loop controller: column vector dimension
                        (sys.input_dim*N) x 1
    risk_alloc_state 
                - Risk allocation for the state constraints
    risk_alloc_input
                - Risk allocation for the input constraints
 
  See also SReachPoint.
 
  Notes:
  * See @LtiSystem/getConcatMats for more information about the notation used.
  
  ============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```
