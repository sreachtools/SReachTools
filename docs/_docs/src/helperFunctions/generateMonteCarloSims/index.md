---
layout: docs
title: generateMonteCarloSims.m
---

```
  Generate Monte-Carlo simulations for a Gaussian-perturbed LTI/LTV system
  (controlled or uncontrolled)
  ============================================================================
  
  generateMonteCarloSims produces a required number of trajectories for a
  Gaussian LTI system.
 
  See also examples/forwardStochasticReachCWH.m
 
  =============================================================================
  concat_state_realization = generateMonteCarloSims(n_monte_carlo_sims, ...
                                                    sys, ...
                                                    initial_state, ...
                                                    time_horizon, ...
                                                    optimal_input_vector)
 
  Inputs:
  -------
    n_monte_carlo_sims   - Number of Monte-Carlo simulation particles to be used 
                           for estimation of the reach-avoid probability
    sys                  - System description as a LtiSystem/LtvSystem object
    initial_state        - Deterministic x_0
    time_horizon         - Time horizon (N) of the stochastic reach-avoid
                           problem
    optimal_input_vector - [Optional] Open-loop controller, a column vector of
                           dimension (sys.input_dim*N) x 1 | Required only if 
                           the system is controlled
 
 
  Outputs:
  --------
    concat_state_realization  - Matrix of concatenate state (row) vectors
                                stacked columnwise. Each row comprises of
                                the state trajectory as [x_1; x_2; ...; x_N]
    concat_disturb_realization- Matrix of concatenate disturbance (row) vectors
                                stacked columnwise. Each row comprises of
                                the state trajectory as [w_0; w_1; ...; w_{N-1}]
 
  Notes:
  ------
  * Assumes IID Gaussian disturbance for the LTI/LTV system. 
  * For controlled system, an open-loop controller NEEDS to be provided. The
    optimal_input_vector should be a ((sys.input_dim) *
    time_horizon)-dimensional vector U = [u_0; u_1; ...; u_N] (column vector).
  * For uncontrolled system, the optimal_input_vector NEED NOT be provided
  
  ============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```
