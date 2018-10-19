---
layout: docs
title: SReachSetLag.m
---

```
  Get approximate level set using lagrangian methods
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
    sys              - LtiSystem object
    prob_thresh             - Probability threshold
    safety_tube      - Cell array of polyhedron objects
    approx_type      - Approximation type, either 'overapproximation' or
                       'underapproximation'
    method, varargin - See help getBoundedSetForDisturbance
 
  Outputs:
  --------
    approx_level_set - Polyhedron object
 
  ============================================================================
  
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
  
```
