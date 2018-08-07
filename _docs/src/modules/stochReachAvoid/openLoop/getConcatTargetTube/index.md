---
layout: docs
title: getConcatTargetTube.m
---

```
  SReachTools/stochasticReachAvoid/getConcatTargetTube: Get concatenated
  target tube
  ============================================================================
 
  getConcatTargetTube computes the concatenated target tube,
  safe_set^{time_horizon -1 } x target_set, a huge polyhedron in the
  (sys.state_dim x time_horizon)-dimensional Euclidean space.
  The output matrices satisfy the relation that the a concatenated state vector
  X lies in the reach-avoid tube if and only if
  
  concat_target_tube_A * X <= concat_target_tube_b 
 
  Usage: See getLowerBoundStochReachAvoid.
 
  ============================================================================
 
  [concat_target_tube_A, concat_target_tube_b] = ..
      getConcatTargetTube(safe_set, target_set, time_horizon);
  
  Inputs:
  -------
    time_horizon         - Time horizon of the stochastic reach-avoid problem
    safe_set             - Safe set for the stochastic reach-avoid problem
    target_set           - Target set for the stochastic reach-avoid problem
 
  Outputs:
  --------
    concat_target_tube_A - State matrix concatenated for target tube
    concat_target_tube_b - Input matrix concatenated for target tube
 
  Notes:
  ------
  * This function also serves as a delegatee for input handling.
  
  ============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
 
```
