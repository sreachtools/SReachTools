---
layout: docs
title: SReachPointGpO.m
---

```
  Solve the stochastic reach-avoid problem (lower bound on the probability and 
  an open-loop controller synthesis) using Genz's algorithm and MATLAB's
  patternsearch (a nonlinear, derivative-free, constrained optimization solver)
  =============================================================================
 
  SReachPointGpO implements the Fourier transform-based underapproximation to
  the stochastic reachability of the target tube problem. A simpler reach-avoid
  problem formulation was discussed in
 
  A. Vinod and M. Oishi, "Scalable Underapproximation for Stochastic
  Reach-Avoid Problem for High-Dimensional LTI Systems using Fourier
  Transforms," in IEEE Control Systems Letters (L-CSS), 2017.
 
  =============================================================================
 
   [lb_stoch_reach, opt_input_vec] = SReachPointGpO(sys, initial_state, ...
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
    lb_stoch_reach 
                - Lower bound on the stochastic reachability of a target tube
                    problem computed using convex chance constraints and
                    piecewise affine approximation
    opt_input_vec
                - Open-loop controller: column vector of dimension
                  (sys.input_dim*N) x 1
 
  See also SReachPoint.
 
  Notes:
  ------
  * Uses SReachPoint('term','chance-open', ...) for initialization of
    patternsearch
  * Uses Genz's algorithm (see in src/helperFunctions) instead of MATLAB's
    Statistics and Machine Learning Toolbox's mvncdf to compute the integral of
    the Gaussian over a polytope
  * See @LtiSystem/getConcatMats for more information about the notation used.
  
  ============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```
