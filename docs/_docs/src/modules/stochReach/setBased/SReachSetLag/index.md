---
layout: docs
title: SReachSetLag.m
---

```
  Get approximate stochastic reach set using Lagrangian methods
  ============================================================================
 
  This function will get the approximate prob_thresh level set for a stochastic
  discrete time system using the Lagrangian methods in 
      J. D. Gleason, A. P. Vinod, M. M. K. Oishi, "Underapproximation of 
      Reach-Avoid Sets for Discrete-Time Stochastic Systems via Lagrangian 
      Methods," in Proceedings of the IEEE Conference on Decision and Control, 
      2017
 
  Usage: see examples/doubleIntegratorLevelSetApprox.m
         or  examples/lagrangianApproximations.m
  
  ============================================================================
 
  Inputs:
  -------
    method_str  - Lagrangian method,
                    'lag-over'  -- Lagrangian Overapproximation
                    'lag-under' -- Lagrangian Underapproximation
    sys         - LtiSystem object
    prob_thresh - Probability threshold
    safety_tube - Tube object
    options     - Struct of reach set options, see SReachSetOptions
 
  Outputs:
  --------
    approx_set - Polyhedron object for the over-/under-approximation of the 
                 stochastic reach set
    approx_tube- [Optional] Tube comprising of an over-/under-approximation of
                 the stochastic reach sets across the time horizon
    bounded_set- [Optional] Bounded disturbance set which was used to
                 robustify the computation (lag-under) or 
                 augment the input set (lag-over)
 
  Notes:
  ------
  * compute_style of `support` and method of `lag-over` will return only the 
    approx_set.
  * While 'Gaussian' disturbance can have options.bound_set_method be 'polytope'
    or 'ellipsoid', 'UserDefined' disturbance requires options.bound_set_method
    to be 'polytope'.
  ============================================================================
  
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
  
```
