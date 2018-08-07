---
layout: docs
title: computeFtLowerBoundStochReachAvoid.m
---

```
  SReachTools/stochasticReachAvoid/computeFtLowerBoundStochReachAvoid: Solve 
  the stochastic reach-avoid problem (lower bound on the probability and an 
  open-loop controller synthesis) using Fourier transform and convex 
  optimization (Internal function --- assumes arguments are all ok)
  =============================================================================
 
  computeFtLowerBoundStochReachAvoid implements the Fourier
  transform-based underapproximation to the terminal hitting-time stochastic
  reach-avoid problem discussed in
 
  A. Vinod and M. Oishi, "Scalable Underapproximation for Stochastic
  Reach-Avoid Problem for High-Dimensional LTI Systems using Fourier
  Transforms," in IEEE Control Systems Letters (L-CSS), 2017.
 
  USAGE: This function is intended for internal use as it does not sanitize the
  inputs. Please use getLowerBoundStochReachAvoid instead.  For the use of
  user-provided initial guess, see getUnderapproxStochReachAvoidSet.
 
  =============================================================================
  [lower_bound_stoch_reach_avoid, optimal_input_vector] = ...
               computeFtLowerBoundStochReachAvoid(sys, ...
                                                  time_horizon, ...
                                                  concat_input_space_A, ... 
                                                  concat_input_space_b, ...
                                                  concat_target_tube_A, ... 
                                                  concat_target_tube_b, ...
                                                  H, ...
                                                  mean_X_sans_input, ...
                                                  cov_X_sans_input, ...
                                                  guess_optimal_input_vector, ...
                                                  desired_accuracy, ...
                                                  PSoptions)
  
  Inputs:
  -------
    sys                         - LtiSystem object describing the system to be
                                  verified
    time_horizon                - Time horizon of the stochastic reach-avoid
                                  problem
    concat_input_space_A,       
     concat_input_space_b       - (A,b) Halfspace representation for the
                                   polytope U^{time_horizon} set.        
    concat_target_tube_A,       
     concat_target_tube_b       - (A,b) Halfspace representation for the
                                  target tube. For example, the terminal
                                  reach-avoid problem requires a polytope of the
                                  form safe_set^{time_horizon-1} x target_set.        
    H                           - Concatenated input matrix (see
                                  @LtiSystem/getConcatMats for the
                                  notation used)
    mean_X_sans_input           - Mean of X
    cov_X_sans_input            - Covariance of X
    guess_optimal_input_vector  - User provided initial guess for optimal input
                                  vector [Use '[]' if unavailable]
    desired_accuracy            - Accuracy expected for the integral of the
                                  Gaussian random vector X over the
                                  concatenated_target_tube [Use 5e-3 if unsure]
    PSoptions                   - Options for patternsearch [Use '[]' if unsure]
 
  Outputs:
  --------
    lower_bound_stoch_reach_avoid - Lower bound on the terminal-hitting 
                                    stochastic reach avoid problem computed 
                                    using Fourier transform and convex 
                                    optimization
    optimal_input_vector          - Optimal open-loop policy
                                    ((sys.input_dim) *
                                    time_horizon)-dimensional vector 
                                    U = [u_0; u_1; ...; u_N] (column vector)
 
  See also verificationOfCwhDynamics, getLowerBoundStochReachAvoid,
  getUnderapproxStochReachAvoidSet, computeReachAvoidProb
 
  Notes:
  ------
  * NOT ACTIVELY TESTED: Builds on other tested functions.
  * MATLAB DEPENDENCY: Uses MATLAB's Global Optimization Toolbox; Statistics and
                       Machine Learning Toolbox.
                       Needs patternsearch for gradient-free optimization
                       Needs normpdf, normcdf, norminv for Genz's algorithm
  * EXTERNAL DEPENDENCY: Uses CVX (optional)
                       Needs CVX to setup a convex optimization problem that
                       initializes the patternsearch-based optimization. If CVX
                       is unavailable, the user may provide a guess for the
                       initialization.
  * Uses Genz's algorithm (see in src/helperFunctions) instead of MATLAB's
    Statistics and Machine Learning Toolbox's mvncdf to compute the integral of
    the Gaussian over a polytope
  * See @LtiSystem/getConcatMats for more information about the
    notation used.
  
  ============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```
