---
layout: docs
title: polytopesFromMonteCarloSims.m
---

```
  Get/plot convex hulls corresponding to a Monte-Carlo simulation
  ============================================================================
  
  Computes (and plots) the convex hull of the Monte-Carlo simulated trajectories 
  at each time step within the horizon. This function is typically useful in 
  conjunction with generateMonteCarloSims to understand the spread of the 
  Monte-Carlo simulated trajectories.
 
  Usage: See examples/dubinsSReachPoint.m
 
  =============================================================================
 
  set_of_ellipsoids = polytopesFromMonteCarloSims(...
     concat_state_realization, state_dim, relv_states, plot_options)
 
  Inputs:
  -------
    concat_state_realization  
                      - Matrix of concatenated state (column) vectors stacked 
                        columnwise. Each column has the state trajectory 
                        [x_1; x_2; ...; x_N]
    state_dim         - State dimension
    relv_states       - A two-element vector with indices of relevant states
                        among the states indexed from 1 to state_dim
    plot_options      - Plot options (MPT3 Polyhedron/plot) that are directly 
                        passed to the plot. Leave empty if plotting is undesired
 
  Outputs:
  --------
    set_of_polytopes  - Set of polytopes comprised of the points. (Array)
 
  See also generateMonteCarloSims
 
  Notes:
  ------
  * Requires MPT3 to compute the convex hull
  * Note that the initial state is NOT a part of the 
    concatenated_state_realization
  ============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```
