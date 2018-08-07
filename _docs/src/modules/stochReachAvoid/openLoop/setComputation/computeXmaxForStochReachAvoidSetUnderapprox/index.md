---
layout: docs
title: computeXmaxForStochReachAvoidSetUnderapprox.m
---

```
  SReachTools/stochasticReachAvoid/computeXmaxForStochReachAvoidSetUnderapprox: 
  Computes the maximum attainable stochastic reach-avoid probability using 
  open-loop controllers (Internal function --- assumes arguments are all ok)
  ====================================================================================
 
  computeXmaxForStochReachAvoidSetUnderapprox computes the initial
  state x_max and its associated open-loop vector that maximizes the Fourier
  transform-based underapproximation to the terminal hitting-time stochastic
  reach-avoid problem discussed in
 
  A. Vinod, and M. Oishi, "Scalable Underapproximative Verification of
  Stochastic LTI Systems Using Convexity and Compactness," in Proceedings of
  Hybrid Systems: Computation and Control (HSCC), 2018. 
 
  It first initializes for an xmax that is deepest with a feasible input_vector
  that keeps the mean trajectory within the reach-avoid tube. This initial guess
  is then refined using patternsearch.
 
  USAGE: This function is not intended for public use. See
  getUnderapproxStochReachAvoidSet
 
  ==============================================================================
  [maximum_underapproximate_reach_avoid_probability, ...
   xmax, ...
   optimal_input_vector_for_xmax] = ...
                     computeXmaxForStochReachAvoidSetUnderapprox(...
                                      sys, ...
                                      time_horizon, ...
                                      init_safe_set, ...
                                      concat_input_space_A, ... 
                                      concat_input_space_b, ...
                                      concat_target_tube_A, ... 
                                      concat_target_tube_b, ...
                                      Abar, ...
                                      H, ...
                                      mean_X_sans_input_sans_initial_state, ...
                                      cov_X_sans_input, ...
                                      affine_hull_of_interest_2D, ...
                                      desired_accuracy, ...
                                      PSoptions)
 
  Inputs:
  -------
   sys                          - LtiSystem object describing the system to be
                                  verified
   time_horizon                 - Time horizon of the stochastic reach-avoid
                                  problem
   initial_state                - Initial state of interest
   init_safe_set                - Safe set for initial state. It should
                                  also include (if any) the affine equality
                                  constraints such that xmax +
                                  set_of_direction vectors span this set.
                                  [Example: Polyhedron('He',[A_eq, b_eq])]%  concat_input_space_A,  
    concat_input_space_b        - (A,b) Halfspace representation for the
                                  polytope U^{time_horizon} set.        
   concat_target_tube_A,  
    concat_target_tube_b        - (A,b) Halfspace representation for the
                                  target tube. For example, the terminal
                                  reach-avoid problem requires a polytope of the
                                  form safe_set^{time_horizon-1} x target_set.        
   Z                            - Concatenated state matrix (see
                                  @LtiSystem/getConcatMats for the
                                  notation used in next three inputs)
   H                            - Concatenated input matrix
   mean_X_sans_input_sans_initial_state
                                - Mean of X with zero input under the
                                  disturbance and with zero initial state
                                  influence
   cov_X_sans_input             - Covariance of X with zero input under the
                                  disturbance
   desired_accuracy             - Accuracy
   PSoptions                    - Options for patternsearch 
 
  Outputs:
  --------
   maximum_underapproximate_reach_avoid_probability
                                - Maximum terminal hitting-time stochastic
                                  reach-avoid probability that may be attained
                                  via an open-loop controller
   xmax                         - Initial state that attains this maximum
   optimal_input_vector_for_xmax- Optimal input vector for xmax
 
  See also verificationOfCwhDynamics, getUnderapproxStochReachAvoidSet,
  computeReachAvoidProb
 
  Notes:
  ------
  * NOT ACTIVELY TESTED: Builds on other tested functions.
  * MATLAB DEPENDENCY: Uses MATLAB's Global Optimization Toolbox; Statistics and
                       Machine Learning Toolbox.
                       Needs patternsearch for gradient-free optimization
                       Needs normpdf, normcdf, norminv for Genz's algorithm
  * EXTERNAL DEPENDENCY: Uses MPT3 and CVX
                       Needs MPT3 for defining a controlled system and the
                       definition of the safe and the target (polytopic) sets
                       Needs CVX to setup convex optimization problems that
                       initializes the patternsearch-based optimization
  * See @LtiSystem/getConcatMats for more information about the
      notation used.
 
  ==============================================================================
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
 
```
