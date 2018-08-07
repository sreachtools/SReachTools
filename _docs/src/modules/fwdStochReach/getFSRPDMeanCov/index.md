---
layout: docs
title: getFSRPDMeanCov.m
---

```
  SReachTools/forwardStochasticReach/getFSRPDMeanCov: Compute the mean and the
  covariance of the state at a time instant in future
  ============================================================================
 
  Computes the mean and the covariance of a Gaussian-perturbed LTI uncontrolled
  system. This function implements Proposition 1 of
 
  A. Vinod, B. HomChaudhuri, and M. Oishi, "Forward Stochastic Reachability
  Analysis for Uncontrolled Linear Systems using Fourier Transforms", In
  Proceedings of the 20th International Conference on Hybrid Systems:
  Computation and Control (HSCC), 2017.
 
  Usage: See getProbReachSet, examples/forwardStochasticReachCWH.mlx.
 
  ============================================================================
  
  [mean_x, cov_x] = getFSRPDMeanCov(sys, initial_state, target_time)
  
  Inputs:
  -------
    sys           - An object of LtiSystem class 
    initial_state - Initial state can be a deterministic n-dimensional vector
                    or a RandomVector object
    target_time   - Time of interest (positive scalar)
 
  Outputs:
  --------
    mean_x        - Mean of the stochastic disturbance
    cov_x         - Covariance of the stochastic disturbance
 
  See also getProbReachSet, getHmatMeanCovForXSansInput
 
  Notes:
  ------
  * getHmatMeanCovForXSansInput computes the FSRPD for the joint state vector
    and getFSRPDMeanCov computes the FSRPD for the state at time t.
 
  ============================================================================
 
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
 
```
