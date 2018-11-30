---
layout: docs
title: getSReachLagOverapprox.m
---

```
  Get the overapproximation of the stoch reach set
  ============================================================================
 
  This function will compute the overapproximation of the stochastic reach
  set via Algorithm 2 in
  
       J. D. Gleason, A. P. Vinod, and M. M. K. Oishi. 2018. Lagrangian 
       Approximations for Stochastic Reachability of a Target Tube. 
       online. (2018). https://arxiv.org/abs/1810.07118
 
  Usage: See examples/lagrangianApproximations.m
    
  ============================================================================
 
  [overapprox_set, overapprox_tube] = getSReachLagUnderapprox(sys,...
        target_tube, disturbance_set)
 
  Inputs:
  -------
    sys             - LtiSystem object
    target_tube     - Tube object 
    disturbance_set - Polyhedron object (bounded disturbance set)
 
  Outputs:
  --------
    overapprox_set - Polyhedron object for the overapproximation of the 
                     stochastic reach set
    overapprox_tube- [Optional] Tube comprising of an overapproximation of the
                     stochastic reach sets across the time horizon
 
  Notes:
  * From computational geometry, intersections and Minkowski differences are
    best performed in facet representation and Minkowski sums are best
    performed in vertex representation. However, since in this computation,
    all three operations are required, scalability of the algorithm is severly
    hampered, despite theoretical elgance.
 
  ============================================================================
 
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
```
