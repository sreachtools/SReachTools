---
layout: docs
title: getAugEffTarget.m
---

```
  Get augmented effective target set
  ============================================================================
 
  This function will compute the augmented effect target via the algorithm in
  the paper:
       [[Will fill out this once paper is actually submitted]]
 
  Usage: See examples/lagrangianApproximations.m
    
  ============================================================================
 
  Inputs:
  -------
    sys          - LtiSystem object
    target_tube  - Cell array of Polyhedron objects 
    disturbance  - Polyhedron object (bounded disturbance set)
 
  Outputs:
  --------
    aug_eff_target - Polyhedron object for the augmented effective
                                 target set
 
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
