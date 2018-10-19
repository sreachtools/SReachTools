---
layout: docs
title: generateMonteCarloSims.m
---

```
  Generate Monte-Carlo simulations for a Gaussian LTI system (controlled or 
  uncontrolled)
  ============================================================================
  
  generateMonteCarloSims produces a required number of trajectories for a
  Gaussian LTI system.
 
  Usage: See checkViaMonteCarloSims and examples/forwardStochasticReachCWH.mlx
 
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
    sys                  - LtiSystem object describing the system to be verified
    initial_state        - Deterministic x_0
    time_horizon         - Time horizon (N) of the stochastic reach-avoid
                           problem
    optimal_input_vector - (Optional) Optimal open-loop policy. Required only if 
                           the system is controlled
 
 
  Outputs:
  --------
    concat_state_realization  - Matrix of concatenate state (row) vectors
                                stacked columnwise. Each row comprises of
                                the state trajectory as [x_1; x_2; ...; x_N]
    concat_disturb_realization- Matrix of concatenate disturbance (row) vectors
                                stacked columnwise. Each row comprises of
                                the state trajectory as [w_0; w_1; ...; w_{N-1}]
 
  See also checkViaMonteCarloSims
 
  Notes:
  ------
  * MATLAB DEPENDENCY: Uses MATLAB's Statistics and Machine Learning Toolbox
    (mvnrnd)
  * INPUT HANDLING: Delegates part of input handling to @LtiSystem/getConcatMats
  * Assumes IID Gaussian disturbance for the LTI system. 
  * For uncontrolled system, the optimal_input_vector NEED NOT be provided
  * For controlled system, an open-loop controller NEEDS to be provided. The
    optimal_input_vector should be a ((sys.input_dim) *
    time_horizon)-dimensional vector U = [u_0; u_1; ...; u_N] (column vector).
  
  ============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```
