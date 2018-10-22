---
layout: docs
title: getSReachLagUnderapprox.m
---

```
  Get underapproximation of stochastic reach set
  =========================================================================
 
  This function will compute the underapproximation of the stochastic reach
  set via Algorithm 1 in
  
       J. D. Gleason, A. P. Vinod, and M. M. K. Oishi. 2018. Lagrangian 
       Approximations for Stochastic Reachability of a Target Tube. 
       online. (2018). https://arxiv.org/abs/1810.07118
 
  Usage: See examples/lagrangianApproximations.m
 
  =========================================================================
 
  underapprox_set = getRobustEffTarget(sys, ...
                                       target_tube, ...
                                       disturbance, ...
                                       Name, Value)
  Inputs:
  -------
    sys          - LtiSystem object
    target_tube  - Tube object
    disturbance  - Polyhedron object (bounded disturbance set)
 
  Outputs:
  --------
    underapprox_set - Polyhedron object
 
  Notes:
  * From computational geometry, intersections and Minkowski differences are
    best performed in facet representation and Minkowski sums are best
    performed in vertex representation. However, since in this computation,
    all three operations are required, scalability of the algorithm is severly
    hampered, despite theoretical elegance.
    
  =========================================================================
  
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
  
```
