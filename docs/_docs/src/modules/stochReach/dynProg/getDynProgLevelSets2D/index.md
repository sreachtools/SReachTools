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
 
  Usage: See example doubleIntegratorDynmaicProgramming.m
 
  ============================================================================
 
  poly_array = getDynProgLevelSets2D(sys, prob_x, prob_lvls)
  
  Inputs:
  -------
    sys         - LtiSystem object (Needs to be 2-dimensional)
    prob_x      - Mx1 Array of probability values at each grid point in
                  grid_X | Use getDynProgSolForTube to compute this vector
    prob_lvls   - A vector containing safety probability thresholds of interest 
                  Each element needs to be within [0,1].
    target_tube - Target tube used for the dynamic programming solution
 
  Outputs:
  --------
 
  Notes:
  ------
  * Uses MATLAB's contour function to create the boundary info which is then fed
    to MPT (takes a convex hull). Because contour function may ignore corners
    (if value function saturates), we check for all the corners when
    constructing the polytope.
  * To be used in conjunction with getDynProgSolForTube
  
  ============================================================================
  
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
```
