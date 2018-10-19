---
layout: docs
title: computeReachProb.m
---

```
  Compute reach prob using Genz's algorithm
  =============================================================================
 
  computeReachProb computes the objective function of the Fourier
  transform-based underapproximation of the terminal hitting-time stochastic
  reach avoid problem as discussed in
 
  A. Vinod and M. Oishi, "Scalable Underapproximation for Stochastic
  Reach-Avoid Problem for High-Dimensional LTI Systems using Fourier
  Transforms," in IEEE Control Systems Letters (L-CSS), 2017.
 
  Specifically, computeReachProb computes the integral of the Gaussian
  random vector (concatenated state vector) X over the reach-avoid (polytopic)
  tube safe_set^{time_horizon-1} x target_set. 
 
  USAGE: See computeFtLowerBoundStochReachAvoid.
 
  =============================================================================
 
  [reach_avoid_prob] = computeReachProb(input_vector, ...
                                             mean_X_sans_input, ...
                                             cov_X_sans_input, ...
                                             H, ...
                                             concat_target_tube_A, ...
                                             concat_target_tube_b, ...
                                             desired_accuracy)
  
  Inputs:
  -------
    input_vector         - Concatenated input vector under investigation
    mean_X_sans_input    - Mean of (X - H * input_vector)
    cov_X_sans_input     - Covariance matrix of X (Since addition of a
                           constant to a Gaussian doesn't affect the
                           covariance matrix)
    H                    - Concatenated H matrix, See 
                           @LtiSystem/getConcatMats.m
    concat_target_tube_A - concatenated target tube polyhedral definition
    concat_target_tube_b - concatenated target tube polyhedral definition
    desired_accuracy     - Accuracy expected for the integral of the
                           Gaussian random vector X over the
                           concatenated_target_tube
 
  Outputs:
  --------
    reach_avoid_prob - Reach-avoid probability attained using the given
                             input_vector
 
  See also iteratedQscmvnv.
 
  Notes:
  ------
  * NOT ACTIVELY TESTED: TODO
  * NO INPUT HANDLING: For computational speed. To be used via
    computeFtLowerBoundStochReachAvoid
  * MATLAB DEPENDENCY: Uses MATLAB's Statistics and Machine Learning Toolbox.
                       Need normpdf, norminv, normcdf for Genz's algorithm
  * Uses Genz's algorithm in an interative manner to compute the integral of a
    Gaussian over a polytope to desired_accuracy provided. See
    helperFunctions/iteratedQscmvnv.m for more details.
  * In the event, the integral is below the desired_accuracy,
    reach_avoid_prob is set to desired_accuracy. This is to allow to take
    log of the reach_avoid_prob.
  
  =============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
 
```
