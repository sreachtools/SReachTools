---
layout: docs
title: SReachPointGpO.m
---

```
  Solve the problem of stochastic reachability of a target tube (a lower bound
  on the maximal reach probability and an open-loop controller synthesis) using
  Genz's algorithm and MATLAB's patternsearch (a nonlinear, derivative-free,
  constrained optimization solver)
  =============================================================================
 
  SReachPointGpO implements the Fourier transform-based underapproximation to
  the stochastic reachability of a target tube problem. The original problem was
  formulated (for the simpler problem of terminal hitting-time stochastic
  reach-avoid problem) in
 
  A. Vinod and M. Oishi, "Scalable Underapproximation for Stochastic
  Reach-Avoid Problem for High-Dimensional LTI Systems using Fourier
  Transforms," in IEEE Control Systems Letters (L-CSS), 2017.
 
     High-level desc.   : Maximize the multivariate Gaussian integral over a
                          polytope, evaluated using Genz's algorithm, and
                          optimize the nonlinear (log-concave) problem using
                          MATLAB's patternsearch
     Approximation      : Approximate upto a user-specified tolerance
     Controller type    : Open-loop controller that satisfies the hard input
                          bounds
     Optimality         : Optimal open-loop controller for the
                          underapproximation problem due to convexity guarantees
 
  =============================================================================
 
   [approx_stoch_reach, opt_input_vec] = SReachPointGpO(sys, initial_state, ...
       safety_tube, options)
 
  Inputs:
  -------
    sys          - System description (LtvSystem/LtiSystem object)
    initial_state- Initial state for which the maximal reach probability must be
                   evaluated (A numeric vector of dimension sys.state_dim)
    safety_tube  - Collection of (potentially time-varying) safe sets that
                   define the safe states (Tube object)
    options      - Collection of user-specified options for 'genzps-open'
                   (Matlab struct created using SReachPointOptions)
 
  Outputs:
  --------
    approx_stoch_reach 
                - Approximation of the stochastic reachability of a target tube
                  problem | Returns -1 if patternsearch fails (exitflag < 1)
    opt_input_vec
                - Open-loop controller: column vector of dimension
                  (sys.input_dim*N) x 1
 
  See also SReachPoint.
 
  Notes:
  ------
  * We recommend using this function through SReachPoint.
  * This function requires MATLAB's Global Optimization Toolbox for its
    nonlinear solver 'patternsearch'.
  * This function requires CVX to work, since it uses SReachPointCcO for
    initialization of MATLAB's 'patternsearch'.
  * This function uses Genz's algorithm (see in src/helperFunctions) instead of
    MATLAB's Statistics and Machine Learning Toolbox's mvncdf to compute the
    integral of the Gaussian over a polytope.
  * See @LtiSystem/getConcatMats for more information about the notation used.
  * This code is also used internally by SReachSetGpO (stochastic reach set
    underapproximation via genzps-open method). There are a few adjustments done 
    for computational purposes:
    1. Using options.thresh, an inner min operation is used that is
       convexity-preserving. This is an attempt to ensure that the quasi
       Monte-Carlo simulation-driven optimization:
       - does report a safety probability that is above prob_thresh, and
       - does not spend too much time looking for global optimality, when a 
         certificate of exceeding a lower bound suffices. 
    2. After the optimization, the optimal value is reevaluated using a fresh
       set of particles for generality.
  
  ============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```
