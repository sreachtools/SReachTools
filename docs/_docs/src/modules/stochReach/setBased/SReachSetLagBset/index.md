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
  
  bounded_set = getBoundedSetForDisturbance(disturbance, time_horizon,...
    onestep_prob_thresh, option)
                                         e.g.
            bounded_set = getBoundedSetForDisturbance(...
                RandomVector('Gaussian', zeros(2,1), eye(2)), ...
                4, ...
                0.8, ...
                'random', ...
                100);
            bounded_set = getBoundedSetForDisturbance(...
                RandomVector('Gaussian', zeros(2,1), eye(2)), ...
                4, ...
                0.8, ...
                'box');
            bounded_set = getBoundedSetForDisturbance(...
                [], ...
                [], ...
                [], ...
                'load', ...
                '/path/to/the/file/to/load/file.mat');
 
  
  Inputs:
  -------
    disturbance - RandomVector object
    time_horizon- Length of the time horizon
    onestep_prob_thresh - Probability threshold
    option      - TODO
 
  Outputs:
  --------
    bounded_set    - Polyhedron object
 
  Notes:
    - When using the 'load' method the mat files must have only one variable
      saved in the mat file and that variable must be a Polyhedron object.
 
  ============================================================================
 
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
```
