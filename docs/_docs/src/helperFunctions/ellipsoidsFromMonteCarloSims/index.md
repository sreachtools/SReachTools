---
layout: docs
title: ellipsoidsFromMonteCarloSims.m
---

```
  Get/plot (minimum volume) ellispoids corresponding to a Monte-Carlo simulation
  ============================================================================
  
  Computes (and plots) the Lowner-John ellipsoid (Boyd and Vandenberghe, Convex
  Optimization, Ch. 8.4.1) that fits the set of samples at each time instant in
  the trajectory. This function is typically useful in conjunction with
  generateMonteCarloSims to understand the spread of the Monte-Carlo simulated
  trajectories.
 
  Usage: See examples/cwhSReachPoint.m
 
  =============================================================================
 
  set_of_ellipsoids = ellipsoidsFromMonteCarloSims(...
     concat_state_realization, state_dim, relv_states, plot_options)
 
  Inputs:
  -------
    concat_state_realization  
                      - Matrix of concatenate state (row) vectors stacked
                        columnwise. Each row comprises of the state trajectory
                        as [x_1; x_2; ...; x_N]
    state_dim         - State dimension
    relv_states       - A two-element vector with indices of relevant states
                        among the states indexed from 1 to state_dim
    plot_options      - Plot options that are directly passed to the plot. Leave
                        empty if plotting is not desired.
 
  Outputs:
  --------
    set_of_ellipsoids - Set of ellipsoids that cover these points. A cell array
                        of cells are provided with each cell containing the
                        ellipsoid center q and its shape matrix Q. The ellipsoid
                        is (x-q)^T Q^{-1} (x-q) <= 1.
 
  See also generateMonteCarloSims
 
  Notes:
  ------
  * Requires CVX for solving the convex optimization problem using SDPT3
  * Note that the initial state is NOT a part of the 
    concatenated_state_realization
  * Uses code obtained from
    http://web.cvxr.com/cvx/examples/cvxbook/Ch08_geometric_probs/min_vol_elp_finite_set.m
  
  ============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```
