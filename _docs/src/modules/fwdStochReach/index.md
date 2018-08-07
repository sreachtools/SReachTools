---
layout: docs
title: getProbReachTargetTube.m
---

```
  SReachTools/forwardStochasticReach/getProbReachTargetTube: Compute the
  probability that the state will lie in a target tube. The starting point may
  be a vector or a RandomVector object
  ============================================================================
 
  This function uses getHmatMeanCovForXSansInput to compute the forward
  stochastic reach probability density (FSRPD) of the concatenated state vector.
  Next, it evaluates the integral of the resulting Gaussian over the
  user-specified target tube (a MPT Polyhedron) using iteratedQscmvnv.
 
  Usage: See examples/forwardStochasticReachCWH.mlx
 
  ============================================================================
  
  prob = getProbReachTargetTube(sys, ...
                                initial_state, ...
                                target_tube, ...
                                desired_accuracy, ...
                                varargin)
 
  Inputs:
  -------
    sys              - An object of LtiSystem class 
    initial_state    - Initial state can be a deterministic n-dimensional vector
                       or a RandomVector object
    target_tube      - Target tube to stay within [TargetTube object]
    desired_accuracy - Accuracy of the integral evaluation 
                       [Default 1e-3 otherwise]
    input_policy     - (Required only for controlled systems) Input policy 
 
  Outputs:
  --------
    prob             - Probability that the system accomplishes the reach-avoid
                       objective
 
  See also iteratedQscmvnv, getFSRPDMeanCov.
 
  Notes:
  ------
  * In case, the target set is a hyper-cuboid and the state_dim < 25,
    then use mvncdf instead.
  * The safe set and the target set must be Polyhedron objects.
  ============================================================================
 
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
 
```
