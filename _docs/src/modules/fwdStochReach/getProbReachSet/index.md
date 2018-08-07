---
layout: docs
title: getProbReachSet.m
---

```
  SReachTools/forwardStochasticReach/getProbReachSet: Compute the probability
  that the state at target_time will lie in the target_set. The starting point
  may be a vector or a RandomVector object
  ============================================================================
 
  This function uses getFSRPDMeanCov to compute the forward stochastic
  reach probability density (FSRPD) at a future target time. Next, it evaluates
  the integral of the resulting Gaussian over the user-specified target_set (a
  MPT Polyhedron) using iteratedQscmvnv.
 
  Usage: See examples/forwardStochasticReachCWH.mlx
 
  ============================================================================
  
  prob = getProbReachSet(sys, ...
                         initial_state, ...
                         target_set, ...
                         target_time, ...
                         desired_accuracy)
 
  Inputs:
  -------
    sys              - An object of LtiSystem class 
    initial_state    - Initial state can be a deterministic n-dimensional vector
                       or a RandomVector object
    target_time      - Time of interest (positive scalar)
    target_set       - Target set 
    desired_accuracy - Accuracy of the integral evaluation 
                       [If unsure, use 1e-8 if sys.state_dim <= 4
                                       1e-3 otherwise]
 
  Outputs:
  --------
    prob             - Probability that x_{target_time} lies in target_set
 
  See also iteratedQscmvnv, getFSRPDMeanCov.
 
  Notes:
  ------
  * In case, the target set is a hyper-cuboid and the state_dim < 25,
    then use mvncdf instead.
  * The target set must be a Polyhedron object.
  ============================================================================
 
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
 
```