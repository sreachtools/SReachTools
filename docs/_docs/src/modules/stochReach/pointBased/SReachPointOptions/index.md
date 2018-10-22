---
layout: docs
title: SReachPointOptions.m
---

```
  Create user-specifiable options for use with SReachPoint()
  =============================================================================
 
  SReachPointOptions creates a MATLAB struct that contains user-specifiable
  options that may be used with SReachPoint
 
  =============================================================================
 
    options = SReachPointOptions(prob_str, method_str, varargin)
  
  Inputs:
  -------
    prob_str    - String specifying the problem of interest. For each case, we
                  compute the optimal value function that maps initial states
                  to different maximal reach probabilities
                      1. 'term' : Stay within the safety_tube
    method_str  - Solution technique to be used (user-specifiable
                  options associated with each technique is enumerated)
                      'chance-open'  -- Convex chance-constrained approach for
                                        an open-loop controller synthesis
                                        1. pwa_accuracy: 
                                                Accuracy of the piecewise affine 
                                                approximation of norminvcdf
                                                used
                      'chance-affine'-- Convex chance-constrained approach for
                                        an affine controller synthesis
                                        1. verbose: Verbosity of the 
                                                implementation (feedback for the 
                                                user)
                                        2. pwa_accuracy: 
                                                Accuracy of the piecewise affine 
                                                approximation of norminvcdf
                                                used
                                        3. max_input_viol_prob:
                                                Probabilistic relaxation of the 
                                                hard input constraints
                                        Difference-of-convex parameters: 
                                        4. tau_initial: Initialization of the 
                                                slack multiplier
                                        5. scaling_tau: Scaling factor to the 
                                                slack multiplier
                                        6. tau_max: Maximum value for the 
                                                scaling factor
                                        7. iter_max: Maximum number of
                                                iterations for the difference of
                                                convex iterative algorithm
                                        8. dc_conv_tol: Tolerance for exiting 
                                                the iterative algorithm
                                        9. slack_tol: Tolerance for the sum
                                                of slack vars for penalty DC
                      'genzps-open'  -- Genz's algorithm + Patternsearch
                                        1. desired_accuracy: 
                                                Accuracy of Gaussian
                                                integral => Accuracy of the
                                                result
                                        2. PSoptions: 
                                                MATLAB struct generated
                                                using psoptimset()
                      'scenario-open'-- Scenario-based 
 
  Outputs:
  --------
    options     - Collection of user-specified options for 'chance-affine'
                  (Matlab struct created using SReachPointOptions)
 
  See also SReachPoint.
 
  Notes:
  * SReachPoint() will call this function internally using the default
    values if SReachPointOptions()-based options is not explicitly provided
    to SReachPoint().
  ============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```
