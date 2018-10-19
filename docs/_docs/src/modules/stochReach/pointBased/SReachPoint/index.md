---
layout: docs
title: SReachPoint.m
---

```
  Solve the stochastic (first/terminal) reach problem approximately from a given
  initial state using a host of techniques
  =============================================================================
 
  SReachPoint computes an approximation to the first/terminal stochastic reach
  problems. 
 
  This function can (approximately) solve two stochastic reachability
  problems:
 
  1. First hitting-time stochastic reachability problem:
 
      maximize Prob( \cup_{i=1}^N {\cap_{t=0}^{i-1} 
                                x_t lies in Safe_t\TargetHyp, x_i\in\TargetHyp})
      subject to
            dynamics and bounds on control
 
  2. Terminal hitting-time stochastic reachability problem (stochastic
  reachability of a target tube):
 
      maximize Prob( \cap_{i=1}^N x_t lies in Safe_t)
      subject to
            dynamics and bounds on control
 
  Using the theory discussed in,
  
  A. P. Vinod and M. Oishi, HSCC 2019 TODO
 
  We can underapproximate the first hitting-time problem by computing a
  finite maximum of a series of stochastic reachability of a target tube with 
  varying time-horizon. 
 
  For the affine controller synthesis problem, we relax hard bounds on the
  control to a user-specified bound on the probability that the affine
  controller
 
  This function is a compilation of various techniques proposed in the
  literature:
 
  1. Convex chance-constrained-based approach (chance-open):
 
     High-level desc.   : Use Boole's inequality, Gaussian random vector, and
                          piecewise linear approximation of the inverse of the
                          standard normal cumulative density function to create
                          a linear program-based approximation to the original
                          optimization
     Approximation      : Guaranteed underapproximation
     Controller type    : Open-loop controller that satisfies the hard
                          input bounds 
     Optimality         : Optimal open-loop controller for the
                          underapproximation problem due to convexity guarantees
     SReachTool function: SReachPointCcO
     Dependency (EXT)   : CVX
     Dependency (MATLAB): Symbolic toolbox
     Paper              : a. Lesser, Oishi, Erwin TODO.
                          b. A. Vinod and M. Oishi, HSCC 2018 TODO
 
  2. Convex chance-constrained-based approach (chance-affine):
 
     High-level desc.   : Use Boole's inequality, Gaussian random vector,
                          hyperbolic constraints-to-second order cone constraint
                          reformulation, and piecewise linear approximation of
                          the inverse of the standard normal cumulative density
                          function to create a second-order cone program-based
                          difference-of-convex optimization problem
     Controller type    : A history-dependent affine controller that satisfies
                          softened input constraints (controller satisfies the
                          hard input bounds upto a user-specified probabilistic
                          threshold)
     Optimality         : Suboptimal affine controller for the
                          underapproximation problem due to the use of
                          difference-of-convex
     Approximation      : Guaranteed underapproximation
     SReachTool function: SReachPointCcA
     Dependency (EXT)   : CVX
     Dependency (MATLAB): Symbolic toolbox
     Paper              : A. Vinod and M. Oishi, HSCC 2018 TODO
 
  3. Fourier transform + Patternsearch (genzps-open):
 
     High-level desc.   : Maximize the multivariate Gaussian integral over a
                          polytope, evaluated using Genz's algorithm, and
                          optimize the nonlinear (log-concave) problem using
                          MATLAB's patternsearch
     Approximation      : Approximate upto a user-specified tolerance
     Controller type    : Open-loop controller that satisfies the hard input
                          bounds
     Optimality         : Optimal open-loop controller for the
                          underapproximation problem due to convexity guarantees
     Dependency (MATLAB): Global Optimization toolbox (for patternsearch)
     SReachTool function: SReachPointGpO
     Paper              : A. Vinod and M. Oishi, "Scalable Underapproximation
                          for Stochastic Reach-Avoid Problem for
                          High-Dimensional LTI Systems using Fourier
                          Transforms," in IEEE Control Systems Letters, 2017.
 
  4. Scenario-based approach (scenario-open):
 
     High-level desc.   : Sample scenarios based on the additive noise and solve
                          a mixed-integer linear program to make the maximum
                          number of scenarios satisfy the reachability objective
     Approximation      : No direct approximation guarantees. Accuracy improves
                          as the number of scenarios considered increases.
     Controller type    : Open-loop controller that satisfies the hard input
                          bounds
     Optimality         : Optimal (w.r.t scenarios drawn) open-loop controller
                          for the underapproximation problem 
     Dependency (EXT)   : CVX
     SReachTool function: SReachPointScO TODO
     Paper              : Lesser, Oishi, Erwin TODO.
 
 
  USAGE: TODO
 
  =============================================================================
 
  [approx_reach_prob, opt_controller, varargout] = SReachPoint(prob_str,...
     method_str, sys, initial_state, safety_tube, options)
  
  Inputs:
  -------
    prob_str     - String specifying the problem of interest. For each case, we
                   compute the optimal value function that maps initial states
                   to different maximal reach probabilities
                       1. 'first' : Stay within the safety_tube and reach the
                                    target set early if possible
                       2. 'term' : Stay within the safety_tube
    method_str   - Solution technique to be used.
                       'chance-open'  -- Convex chance-constrained approach for
                                         an open-loop controller synthesis
                       'chance-affine'-- Convex chance-constrained approach for
                                         an affine controller synthesis
                       'genzps-open'  -- Genz's algorithm + Patternsearch
                       'scenario-open'-- Scenario-based 
    sys          - System description (LtvSystem/LtiSystem object)
    initial_state- Initial state for which the maximal reach probability must be
                   evaluated (A numeric vector of dimension sys.state_dim)
    safety_tube  - Collection of (potentially time-varying) safe sets that
                   define the safe states (Tube object)
    options      - Collection of user-specified options for each of the solution
                   (Matlab struct created using SReachPointOptions)
 
  Outputs:
  --------
    approx_reach_prob 
                - Approximation (underapproximation, in some cases) to the
                  first/terminal stochastic reach problem
    opt_input_vec, 
      opt_input_gain
                - Controller U=MW+d for a concatenated input vector 
                    U = [u_0; u_1; ...; u_{N-1}] and concatenated disturbance
                    vector W=[w_0; w_1; ...; w_{N-1}]. 
                    - opt_input_gain: Affine controller gain matrix of dimension
                        (sys.input_dim*N) x (sys.dist_dim*N)
                    - opt_input_vec: Open-loop controller: column vector dimension
                        (sys.input_dim*N) x 1
                  The feedback gain matrix M is set to [] for methods that look
                  for open-loop controllers. 
    risk_alloc_state 
                - [Available only for 'chance-X'] Risk allocation for the
                  state constraints
    risk_alloc_input
                - [Available only for 'chance-affine'] Risk allocation for the
                  input constraints
 
  Notes:
  * SReachPoint() will call this function internally using the default
      values if SReachPointOptions()-based options is not explicitly provided
      to SReachPoint().
  * See @LtiSystem/getConcatMats for more information about the notation used.
  * If an open_loop policy is desired arranged in increasing time columnwise,
    use the following command:
        optimal_open_loop_control_policy = reshape(opt_controller,...
            sys.input_dim, time_horizon);
  
  =============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
 
```
