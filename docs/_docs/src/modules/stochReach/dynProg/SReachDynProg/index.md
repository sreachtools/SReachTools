---
layout: docs
title: SReachDynProg.m
---

```
  Dynamic programming solution to stochastic reachability problems
  ============================================================================
 
  The function computes the probability of staying in a target tube defined
  on a particular state stace grid. The dynamic programming recursion can be
  found in 
    
  S. Summers and J. Lygeros, "Verification of discrete time stochastic hybrid 
  systems: A stochastic reach-avoid decision problem," Automatica, vol. 46,
  no. 12, pp. 1951--1961, 2010.
 
  A trivial extension of this work to the case of time-varying safe set is
  implemented here.
 
  Usage: See example/doubleIntegratorDynamicProgramming.mlx
 
  ============================================================================
 
  prob_x = SReachDynProg('term',sys,x_inc,u_inc,safety_tube)
  
  Inputs:
  -------
    prob_str    - String specifying the problem of interest. For each case, we
                  compute the optimal value function that maps initial states
                  to different maximal reach probabilities
                      1. 'term' : Stay within the safety_tube
    sys         - System description as a LtiSystem object
    x_inc       - Scalar increment for all dimensions of the state space
    u_inc       - Scalar increment for all dimensions of the input space
    safety_tube - Safety tube of length N+1 where N is the time_horizon. It
                  should have polyhedrons T_0, T_1, ...,T_N.
 
  Outputs:
  --------
    prob_x      - Probability values at each grid point (M number of them) in
                  grid_X (Mx1 array)
    cell_xvec   - [Optional] Gridding along the particular dimension 
                  (sys.state_dim x 1 cell array, with grid info along each
                  dimension)
    grid_x      - [Optional] Collection of grid points (Mx1 array)
    mat_prob_x  - [Optional] M*(N+1) matrix of probability values corresponding
                  to the "unrolled" value functions [V_0, V_1, ... V_N] where N
                  is the time horizon. Note that prob_x = mat_prob_x(1,:)
 
  See also getDynProgLevelSets2D
 
  Notes:
  ------
  * REQUIRES:
    - Input space is an axis-aligned HYPERCUBOID.
    - State space is the smallest axis-aligned HYPERCUBOID that contains all the
      sets in the target-tube
  * We impose uniform gridding across every dimension.
  * WARNING: Dynamic programming suffers from the curse of dimensionality!
    Using fine grids will increase the computation time.
  
  ============================================================================
  
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
```
