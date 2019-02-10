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
 
  [underapprox_set, underapprox_tube] = getSReachLagUnderapprox(sys,...
        target_tube, disturbance_set)
 
  Inputs:
  -------
    sys              - LtiSystem object
    target_tube      - Tube object
    dist_set         - Polyhedron/SReachEllipsoid object (bounded set) OR a
                        collection of these objects which individually satisfy 
                        the probability bound(a convex hull of the individual 
                        results taken posteriori)
    options          - Struct of reach set options, see SReachSetOptions
 
  Outputs:
  --------
    overapprox_set   - Polyhedron object for the underapproximation of the 
                       stochastic reach set
    underapprox_tube - [Optional] Tube comprising of an underapproximation of
                       the stochastic reach sets across the time horizon
 
  Notes:
  * From computational geometry, intersections and Minkowski differences are
    best performed in facet representation and Minkowski sums are best
    performed in vertex representation. However, since in this computation,
    all three operations are required, scalability of the algorithm is severly
    hampered, despite theoretical elegance.
  * Since box and random approaches in SReachSetOptions produce Polyhedron
    objects for disturbance sets, we rely on MPT for all the set operations.
    This means we do have scalability issues mentioned above.
  * For ellipsoid approach in SReachSetOptions, we seek a purely facet-based
    operation and utilize the ray-shooting algorithm to compute a facet-based
    underapproximation of the Minkowski sum step (via vertex-based
    underapproximation, followed by projection, followed by convex hull
    operation)
  * Use spreadPointsOnUnitSphere.m to compute equi_dir_vecs. 
  * equi_dir_vecs is automatically generated as part of SReachSetOptions.
    
  =========================================================================
  
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
  
```
