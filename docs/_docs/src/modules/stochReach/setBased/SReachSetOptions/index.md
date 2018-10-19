---
layout: docs
title: SReachSetOptions.m
---

```
  Create user-specifiable options for use with SReachSet()
  =============================================================================
 
  SReachSetOptions creates a MATLAB struct that contains user-specifiable
  options that may be used with SReachSet
 
  =============================================================================
 
    options = SReachSetOptions(prob_str, method_str, varargin)
  
  Inputs:
  -------
    prob_str    - String specifying the problem of interest. For each case, we
                  compute the optimal value function that maps initial states
                  to different maximal reach probabilities
                      1. 'term' : Stay within the safety_tube
    method_str  - Solution technique to be used (user-specifiable
                  options associated with each technique is enumerated)
                      'chance-open':
                           Convex chance-constrained approach for an open-loop
                           controller synthesis
                           1. init_safe_set_affine: 
                           1. set_of_dir_vecs: 
                           1. verbose: 
                           1. pwa_accuracy: Accuracy of the piecewise affine
                                approximation of norminvcdf used
                      'genzps-open':
                           Genz's algorithm + Patternsearch
                           1. init_safe_set_affine: 
                           1. set_of_dir_vecs: 
                           1. verbose:
                           1. tol_bisect
                           1. desired_accuracy: Accuracy of Gaussian integral =>
                                Accuracy of the result
                           2. PSoptions: MATLAB struct from psoptimset()
                      'lag-over'/'lag-under':
                           Lagrangian-based over- and underapproximation
                           bound_set_method:
                           1a. 'random'- Get an approximation of the ellipsoid
                                         using random direction choices; only
                                         usable for Gaussian-type disturbances;
                                         varargin must be an integer for the
                                         number of random directions to be used;
                                         b. num_dirs --- Number of directions to
                                            sample the ellipsoid for a polytopic
                                            representation
                           2a. 'box'    - Get an n-dimensional cuboid centered at
                                         the disturbance mean that satisfies the
                                         probability threshold
                                         b. err_thresh --- Tolerance for
                                            the bisection algorithm that
                                            identifies the length of the box
                           3a. 'optim-box' 
                                       - Get an n-dimensional cuboid centered at
                                         a user-specified box center that
                                         satisfies the probability threshold
                                         b. box_center --- Center for the box
                           4a. 'load'   - Load a predefined polyhedron bounding
                                         set; primarily used for comparison and
                                         repeatability testing
                                         b. load_str --- Path to the file to
                                            load. All other inputs are
                                            IRRELEVANT for this option.
                                            Mat files to be loaded must
                                            have specific design TODO
 
  Outputs:
  --------
    options     - Collection of user-specified options for 'chance-affine'
                  (Matlab struct created using SReachSetOptions)
 
  See also SReachSet.
 
  Notes:
  * SReachSet() will call this function internally using the default
      values if SReachSetOptions()-based options is not explicitly provided
      to SReachSet().
  ============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```
