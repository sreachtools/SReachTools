---
layout: docs
title: getLowerBoundStochReachAvoid.m
---

```
  SReachTools/stochasticReachAvoid/getLowerBoundStochReachAvoid: Solve the
  stochastic reach-avoid problem (lower bound on the probability and an
  open-loop controller synthesis) using Fourier transform and convex
  optimization
  =============================================================================
 
  getLowerBoundStochReachAvoid computes a lower bound to the terminal
  hitting-time stochastic reach-avoid problem by searching over the space of
  open-loop controllers.
 
  A. Vinod and M. Oishi, "Scalable Underapproximation for Stochastic Reach-Avoid
  Problem for High-Dimensional LTI Systems using Fourier Transforms," in IEEE
  Control Systems Letters (L-CSS), 2017.
 
  This function has input handling and sets up the problem to be solved by
  depending on the user-preference:
    1. computeFtLowerBoundStochReachAvoid, (Option: genzps)
    2. computeCcLowerBoundStochReachAvoidIterRisk, (Option: ccciter)
    3. computeCcLowerBoundStochReachAvoidPwlRisk (Option: cccpwl).
 
  USAGE: TODO
 
  =============================================================================
 
  [lb_stoch_reach_avoid, optimal_input_vector] =...
                                  getLowerBoundStochReachAvoid(sys,...
                                                               initial_state,...
                                                               target_tube,...
                                                               method,...
                                                               varargin)
 
  Inputs:
  -------
    sys                  - LtiSystem object describing the system to be verified
    initial_state        - Initial state of interest
    target_tube          - Target tube to stay within [TargetTube object]
    method               - Method to compute the reach-avoid probability
                              'genzps' -- Genz's algorithm + Patternsearch
                              'cccpwl' -- Piecewise-linear conservative
                                          implementation for convex
                                          chance-constrained reformulation
                              'ccciter'-- Iterative risk allocation approach for
                                          convex chance-constrained
                                          reformulation
                           See Notes for dependencies.
    guess_optimal_input_vector
                         - (Optional) Provide a concatenated guess for the
                           optimal input policy vector in the form of U = [u_0;
                           u_1; ...; u_N]. [If unsure, provide []. This will
                           trigger a CVX-based initialization computation.]
    desired_accuracy     - (Optional) Accuracy  [Default 5e-3]
    PSoptions            - (Optional) Options for patternsearch [Default
                            psoptimset('Display', 'off')]
 
  Outputs:
  --------
    lb_stoch_reach_avoid - Lower bound on the terminal-hitting stochastic reach
                           avoid problem computed using Fourier transform and
                           convex optimization
    optimal_input_vector - Optimal open-loop policy ((sys.input_dim) *
                           time_horizon)-dim.  vector U = [u_0; u_1; ...; u_N]
                           (column vector)
 
  See also computeFtLowerBoundStochReachAvoid, getCcLowerBoundStochReachAvoid.
 
  Notes:
  * NOT ACTIVELY TESTED: Builds on other tested functions.
  * MATLAB DEPENDENCY : Uses MATLAB's Statistics and Machine Learning Toolbox
                        Needs normpdf, normcdf, norminv for Genz's algorithm
  * EXTERNAL DEPENDENCY: Uses MPT3
                         Needs MPT3 for defining a controlled system and the
                         definition of the safe and the target (polytopic) sets
  * Method 'genzps' has the following dependencies
        * MATLAB DEPENDENCY: Uses MATLAB's Global Optimization Toolbox
                             Needs patternsearch for gradient-free optimization
        * EXTERNAL DEPENDENCY: Uses CVX (optional)
                               Needs CVX to setup a convex optimization problem
                               that initializes the patternsearch-based
                               optimization. If CVX is unavailable, the user may
                               provide a guess for the initialization.
        * Specify both desired_accuracy and PSoptions or neither to use the
          defaults 
        * Specify an optional guess_optimal_input_vector to skip the use of CVX
  * Method 'cccpwl' has the following dependencies
        * EXTERNAL DEPENDENCY: Uses CVX 
                               Needs CVX to setup the convex chance-constrained
                               problems
        * Specify a desired_accuracy if required. Else, a default value of 1e-3
            is used.
  * Method 'ccciter' has the following dependencies
        * EXTERNAL DEPENDENCY: Uses CVX 
                               Needs CVX to setup the convex chance-constrained
                               problems
        * Specify a desired_accuracy if required. Else, a default value of 1e-3
            is used.
  * See @LtiSystem/getConcatMats for more information about the
    notation used.
  * If an open_loop policy is desired arranged in increasing time columnwise,
    use the following command
        optimal_open_loop_control_policy = reshape(...
            optimal_input_vector, ...
            sys.input_dim, ...  
            time_horizon);
  
  =============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
 
```
