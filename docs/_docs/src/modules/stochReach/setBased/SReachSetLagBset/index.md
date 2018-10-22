---
layout: docs
title: SReachSetLagBset.m
---

```
  Get bounded disturbance set for approximation
  ============================================================================
 
  This function will get a bounded disturbance set used to compute robust
  reach avoid sets or robust effective target tubes.
 
  Usage: see examples/boundedDisturbanceSets.m
 
  ============================================================================
  
  bounded_set = SReachSetLagBset(disturbance, time_horizon, ...
    onestep_prob_thresh, option)
  
  Inputs:
  -------
    disturbance         - RandomVector object
    time_horizon        - Length of the time horizon
    onestep_prob_thresh - Probability threshold
    option              - Struct specifying methods to obtain the bounded set
                          see SReachSetOptions
 
  Outputs:
  --------
    bounded_set    - Polyhedron object
 
    
    See also SReachSetOptions
  
  Notes:
  * When using the 'load' method the mat files must have only one variable
    saved in the mat file and that variable must be a Polyhedron object.
 
  ============================================================================
 
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
```
