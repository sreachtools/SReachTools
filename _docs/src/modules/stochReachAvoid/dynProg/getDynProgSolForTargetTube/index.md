---
layout: docs
title: getDynProgSolForTargetTube.m
---

```
  SReachTools/stochasticReachAvoid/getDynProgSolForTargetTube Get dynamic 
  programming grid probability for reachability of target tube
  ============================================================================
 
  The function computes the probability of staying in a target tube defined
  on a particular state stace grid. The dynamic programming recursion can be
  found in 
    
  S. Summers and J. Lygeros, "Verification of discrete time stochastic hybrid 
  systems: A stochastic reach-avoid decision problem," Automatica, vol. 46,
  no. 12, pp. 1951--1961, 2010.
 
  The problem of examining the reachability of a target tube can be found in
  a work that we intend to publish soon :D TODO
 
  Usage: See example doubleIntegratorDynmaicProgramming.m
 
  ============================================================================
 
  grid_prob = getDynProgSolForTargetTube(sys, ...
      state_grid, input_grid, target_tube, varargin)
  
  Inputs:
  -------
    sys         - LtiSystem object
    state_grid  - SpaceGrid object
    input_grid  - InputGrid object
    target_tube
                - Target tube of length N+1 where N is the time_horizon. It should have
                    polyhedrons T_0, T_1,...,T_N.
 
  Outputs:
  --------
    grid_prob   - Nx1 Array of probability values, where N is equivalent
                  to size(state_grid, 1)
 
  Notes:
  ------
  * WARNING: Dynamic programming suffers from the curse of dimensionality! As
    such, this code will effective and with reasonable speed compute dynamic
    programming solutions for 2-dimensional systems with Gaussian 
    disturbances. However, for 3-dimensional systems or larger the required
    computation time, and memory, will exponentially grow to the point that 
    the simulation will take longer than it took for me to get my PhD.
  * Currently this back propagation, and subsequently the entire dynamic 
    programming recursion, only works for Gaussian disturbances.
  
  ============================================================================
  
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
  check inputs
```
