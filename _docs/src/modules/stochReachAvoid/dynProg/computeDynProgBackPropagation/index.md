---
layout: docs
title: computeDynProgBackPropagation.m
---

```
  SReachTools/stochasticReachAvoid/computeDynProgBackPropagation Compute the
  dynamic programming back propagation
  ============================================================================
 
  The function computes the one-step back propagation for the dynamic 
  programming recursion. See
  
  S. Summers and J. Lygeros, "Verification of discrete time stochastic hybrid 
  systems: A stochastic reach-avoid decision problem," Automatica, vol. 46,
  no. 12, pp. 1951--1961, 2010.
 
  Usage: See getDynProgSolForTargetTube
  
  ============================================================================
 
  grid_prob = computeDynProgBackPropagation(sys, ...
      state_grid, input_grid, grid_prob, initial_set)
  
  Inputs:
  -------
    sys         - LtiSystem object
    state_grid  - SpaceGrid object
    input_grid  - InputGrid object
    grid_prob   - Nx1 Array of probability values, where N is equivalent
                  to size(state_grid, 1)
    initial_set - Polyhedron object
 
  Outputs:
  --------
    grid_prob - Nx1 Array of probability values, where N is equivalent
                       to size(state_grid, 1)
 
  See also getDynProgSolForTargetTube
 
  Notes:
  ------
  * Currently this back propagation, and subsequently the entire dynamic 
    programming recursion, only works for Gaussian disturbances.
 
  ============================================================================
 
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
```
