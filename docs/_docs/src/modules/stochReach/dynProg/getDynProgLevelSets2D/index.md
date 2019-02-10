---
layout: docs
title: getDynProgLevelSets2D.m
---

```
  SReachTools/stochasticReachAvoid/getDynProgLevelSets2D Get level sets based 
  on the value function returned by getDynProgSolForTube
  ============================================================================
 
  The function computes an array of polytopes based on the results from
  getDynProgSolForTube
 
  See also examples/doubleIntegratorDynamicProgramming.m, SReachDynProg
 
  ============================================================================
 
  poly_array = getDynProgLevelSets2D(sys, prob_x, prob_lvls)
  
  Inputs:
  -------
    cell_xvec  - Gridding along the particular dimension (sys.state_dim x 1
                 cell array, with grid info along each dimension) | Output for
                 SReachDynProg
    prob_x     - Probability values at each grid point (M number of them) in
                 grid_X (Mx1 array)
    prob_lvls  - A vector containing safety probability thresholds of interest
                 Each element needs to be within [0,1].
    safety_tube- Safety tube used for the dynamic programming solution
 
  Outputs:
  --------
    poly_array - Array of Polyhedron objects for each of the level sets
 
  Notes:
  ------
  * Uses MATLAB's contour function to create the boundary info which is then fed
    to MPT (takes a convex hull). Because contour function may ignore corners
    (if value function saturates), we check for all the corners when
    constructing the polytope.
  * To be used in conjunction with SReachDynProg
  
  ============================================================================
  
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
 
```
