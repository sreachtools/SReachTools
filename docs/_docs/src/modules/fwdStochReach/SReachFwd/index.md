---
layout: docs
title: SReachFwd.m
---

```
  Perform forward stochastic reachability analysis of a Gaussian-perturbed
  linear system
  ============================================================================
 
  Perform forward stochastic reachability analysis of a Gaussian-perturbed
  linear system. This function implements ideas from
 
  A. Vinod, B. HomChaudhuri, and M. Oishi, "Forward Stochastic Reachability
  Analysis for Uncontrolled Linear Systems using Fourier Transforms", In
  Proceedings of the 20th International Conference on Hybrid Systems:
  Computation and Control (HSCC), 2017.
 
  See also examples/forwardStochasticReachCWH.m.
 
  ============================================================================
  
  varargout = SReachFwd(prob_str, sys, initial_state, target_time, varargin)
  
  Inputs:
  -------
    prob_str      - String specifying the problem of interest
                        1. 'state-stoch' : Provide mean and covariance at state
                                           at the specified time
                        2. 'state-prob'  : Compute the probability that the
                                           state will lie in a polytope at the
                                           specified time
                        3. 'concat-stoch': Provide mean and covariance of the
                                           concatenated state vector up to a
                                           specified time
                        4. 'concat-prob' : Compute the probability that the
                                           concatenated state vector up to a
                                           specified time lies in the given
                                           target tube
    sys           - System description as a LtiSystem/LtvSystem object
    initial_state - Initial state as a deterministic n-dimensional vector
                    or a RandomVector object
    target_time   - Time of interest (positive scalar)
    target_set/tube  
                  - [Required only for state/concat-prob] Polyhedron/Tube object
                    over which the probability must be computed
    desired_accuracy 
                  - [Required only for state/concat-prob] Accuracy for the
                    integral
    
 
  Outputs:
  --------
    mean_vec      - ['state/concat-stoch'] Mean of the stochastic disturbance
    cov_mat       - ['state/concat-stoch'] Covariance of the stochastic 
                                           disturbance
    prob          - ['state/concat-prob'] Probability of occurence
 
  Notes:
  ------
  * Assumes IID disturbance.
  * The outputs are either (mean_vec, cov_mat) or (prob), depending on the
    method_str
 
  ============================================================================
 
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
 
```
