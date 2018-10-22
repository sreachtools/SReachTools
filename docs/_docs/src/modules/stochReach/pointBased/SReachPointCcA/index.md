---
layout: docs
title: SReachPointCcA.m
---

```
  Solve the problem of stochastic reachability of a target tube (a lower bound
  on the maximal reach probability and an affine controller synthesis) using
  chance-constrained optimization and difference of convex programming
  =============================================================================
 
  SReachPointCcA implements the chance-constrained underapproximation to the
  problem of stochastic reachability of a target tube to construct an affine
  controller. This technique is discussed in detail in the paper,
 
  A. Vinod and M. Oishi. Affine controller synthesis for stochastic reachability
  via difference of convex programming. In Proc. Hybrid Syst.: Comput. & Ctrl.,
  2019. (submitted). https://hscl.unm.edu/affinecontrollersynthesis/
 
     High-level desc.   : Use Boole's inequality, Gaussian random vector,
                          hyperbolic constraints-to-second order cone constraint
                          reformulation, and piecewise linear approximation of
                          the inverse of the standard normal cumulative density
                          function to create a second-order cone program-based
                          difference-of-convex optimization problem
     Controller type    : A history-dependent affine controller that satisfies
                          softened input constraints (controller satisfies the
                          hard input bounds upto a user-specified probabilistic
                          threshold)
     Optimality         : Suboptimal affine controller for the
                          underapproximation problem due to non-convexity
                          established by the difference of convex formulation
     Approximation      : Guaranteed underapproximation
 
  =============================================================================
 
    [lb_stoch_reach, opt_input_vec, opt_input_gain, ...
     risk_alloc_state, risk_alloc_input] = SReachPointCcA(sys,...
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
  * We recommend using this function through SReachPoint.
  * This function requires CVX to work.
  * This function returns a **lower bound to the maximal reach probability under
    hard input constraints**. This lower bound is obtained by a linear
    transformation of the maximal reach probability associated with the
    unsaturated affine controller using the user-specified likelihood threshold
    on the hard input constraints. See Theorem 1 of the paper cited above.
  * See @LtiSystem/getConcatMats for more information about the notation used.
  
  ============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```
