---
layout: docs
title: getBoundedSetForDisturbance.m
---

```
  SReachTools/stochasticReachAvoid/getBoundedSetForDisturbance: Get bounded 
  disturbance set for approximation
  ============================================================================
 
  This function will get a bounded disturbance set used to compute robust
  reach avoid sets or robust effective target tubes.
 
  Usage: see examples/boundedDisturbanceSets.m
 
  ============================================================================
  
  bounded_set = getBoundedSetForDisturbance(disturbance, ...
      horizon_length, beta, method, varargin)
  
  Inputs:
  -------
    disturbance    - StochasticDisturbance object
    horizon_length - Length of the time horizon
    beta           - Probability threshold
    method         - Method for computing bounded set
    varargin       - Dependent upon method chosen, see below
 
    Available methods:
        'random' - Get an approximation of the ellipsoid using random
                   direction choices; only usable for Gaussian-type
                   disturbances; varargin must be an integer for the
                   number of random directions to be used; e.g.
            bounded_set = getBoundedSetForDisturbance(...
                StochasticDisturbance('Gaussian', zeros(2,1), eye(2)), ...
                4, ...
                0.8, ...
                'random', ...
                100);
 
        'box'    - Get an n-dimensional cuboid that satisfies the
                   probability threshold; does not accept varargins;
                   currenlty not implemented; e.g.
            bounded_set = getBoundedSetForDisturbance(...
                StochasticDisturbance('Gaussian', zeros(2,1), eye(2)), ...
                4, ...
                0.8, ...
                'box');
 
        'load'   - Load a predefined polyhedron bounding set; primarily
                   used for comparison and repeatability testing; varargin
                   must be a character array of the path to the file to
                   load; mat files to be loaded must have specific design,
                   see Notes section; when using load method all other 
                   inputs are irrelevant; e.g.
            bounded_set = getBoundedSetForDisturbance(...
                [], ...
                [], ...
                [], ...
                'load', ...
                '/path/to/the/file/to/load/file.mat');
 
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
