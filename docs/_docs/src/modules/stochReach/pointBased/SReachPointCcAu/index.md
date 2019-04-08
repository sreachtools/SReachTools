---
layout: docs
title: SReachPointCcAu.m
---

```
  Solve the problem of stochastic reachability of a target tube (a lower bound
  on the maximal reach probability and an affine controller synthesis) using
  chance-constrained optimization and uniform risk allocation
  =============================================================================
 
  SReachPointCcAuniform implements the chance-constrained underapproximation to
  the problem of stochastic reachability of a target tube to construct an affine
  controller. This technique is inspired from Algorithms 1 and 2 of 
 
  M. Vitus and C. Tomlin, "On feedback design and risk allocation in chance
  constrained control", In Proc. Conf. Dec. & Ctrl., 2011.
 
  In contrast to their original algorithm, we have a chance constraint on
  the input and the state. Further, the lower bound on the reachability (state 
  constraint) depends on how high the input chance constraint satisfaction 
  probability is. Therefore, we perform two levels of bisection
  --- one to maximize the probability of constraint satisfaction for the
  state, and the other to meet the chance constraint on the input. However,
  to save time, we check only for feasibility in the input bisection.
 
  Subsequently, the obtained solution is discounted for input constraint
  violation using Theorem 1 of
 
  A. Vinod and M. Oishi. Affine controller synthesis for stochastic reachability
  via difference of convex programming. In Proc. Conf. Dec. & Ctrl., 2019.
  (submitted). https://hscl.unm.edu/affinecontrollersynthesis/
 
 
  =============================================================================
 
    [lb_stoch_reach, opt_input_vec, opt_input_gain, risk_alloc_state, ...
        risk_alloc_input] = SReachPointCcAu(sys, initial_state, safety_tube, ...
        options)
  
  Inputs:
  -------
    sys          - System description (LtvSystem/LtiSystem object)
    initial_state- Initial state for which the maximal reach probability must be
                   evaluated (A numeric vector of dimension sys.state_dim)
    safety_tube  - Collection of (potentially time-varying) safe sets that
                   define the safe states (Tube object)
    options      - Collection of user-specified options for 'chance-affine-uni'
                   (Matlab struct created using SReachPointOptions)
 
  Outputs:
  --------
    lb_stoch_reach 
                - Lower bound on the stochastic reachability of a target tube
                  problem computed using chance constraints and
                  difference-of-convex techniques
    opt_input_vec, 
      opt_input_gain
                - Controller U=MW+d for a concatenated input vector 
                    U = [u_0; u_1; ...; u_{N-1}] and concatenated disturbance
                    vector W=[w_0; w_1; ...; w_{N-1}]. 
                    - opt_input_gain: Affine controller gain matrix of dimension
                        (sys.input_dim*N) x (sys.dist_dim*N)
                    - opt_input_vec: Open-loop controller: column vector 
                      dimension
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
