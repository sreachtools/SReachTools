---
layout: docs
title: getChainOfIntegLtiSystem.m
---

```
  Get chain of integrators LTI System
  ============================================================================
  
  Get an n-d discrete chain of integrators in the form of an LtiSystem object
  given a discretization time-step T.
 
  Chain of integrators model:
 
    x_{k+1} = A * x_{k} + B * u_{k} + w_{k}
  
    A = [1, T, T^2/2, ... T^n/n!;
         0, 1, T,     ... T^(n-1)/(n-1)!;
         ...
         0, 0, ...    ... T;
         0, 0, ...    ... 1];
   
    B = [T^n/n!;
         T^(n-1)/(n-1)!
         ...
         T];
  
  Usage:
  ------
    % 3-d chain of integrators with U = [-1,1] and no (empty) disturbance
    sys = getChainOfIntegLtiSystem(3, 0.2, ...
        'InputSpace', Polyhedron('lb', -1, 'ub', 1));
 
  ============================================================================
  
  sys = getChainOfIntegLtiSystem(dim, T, input_space, disturb)
  
  Inputs:
  -------
    dim         - Dimensions
    T           - Discretization time step
    input_space - Input space (Polyhedron)
    disturb     - Disturbance object (Polyhedron / RandomVector / empty)
  
  Outputs:
  --------
    sys         - LtiSystem object describing the chain of integrators
  
  =============================================================================
 
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```
