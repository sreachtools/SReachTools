---
layout: docs
title: getDubinsCarLtv.m
---

```
  SReachTools/systems/getDubinsCarLtv: 
  ============================================================================
  
  
  Usage:
  ------
    % 3-d chain of integrators with U = [-1,1] and no (empty) disturbance
    sys = getChainOfIntegLtiSystem(3, 0.2, ...
        Polyhedron('lb', -1, 'ub', 1), ...
        Polyhedron());
 
  ============================================================================
  
  sys = getChainOfIntegLtiSystem(dim, T, input_space, disturb)
  
  Inputs:
  -------
    dim         - Dimensions
    T           - Discretization time step
    input_space - Input space (Polyhedron)
    disturb     - Disturbance object (Polyhedron / StochasticDisturbance)
  
  Outputs:
  --------
    sys - LtiSystem object
  
  =============================================================================
 
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/abyvinod/SReachTools/blob/master/LICENSE
  
 
```
