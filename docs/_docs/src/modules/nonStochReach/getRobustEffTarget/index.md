---
layout: docs
title: getRobustEffTarget.m
---

```
  Get robust Effective Target Set
  =========================================================================
 
  This function will compute the augmented effect target via the algorithm in
  the paper:
       [[Will fill out this once paper is actually submitted]]
 
  TODO
 
  Usage: See examples/lagrangianApproximations.m
 
  =========================================================================
 
  robust_eff_target = getRobustEffTarget(sys, ...
                                         target_tube, ...
                                         disturbance, ...
                                         Name, Value)
  Inputs:
  -------
    sys          - LtiSystem object
    target_tube  - Target tube of length N+1 where N is the time_horizon. It should have
                   polyhedrons T_0, T_1, ...,T_N.
    disturbance  - Polyhedron object (bounded disturbance set)
  
    Name       | Value
    ----------------------------------------
    style      | 'standard', 'vrep'
 
  Outputs:
  --------
    robust_eff_target - Polyhedron object
 
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
