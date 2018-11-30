---
layout: docs
title: generateMonteCarloSims.m
---

```
  Generate Monte-Carlo simulations for (controlled/uncontrolled) LTI/LTV system 
  =============================================================================
  
  generateMonteCarloSims produces a required number of trajectories,
  n_monte_carlo_sims, for a (affine-controlled/uncontrolled) LTI/LTV system sys
  with a deterministic/RandomVector initial_state for a given time_horizon. 
 
  If the system is controlled, then a causal disturbance-feedback affine
  controller may be specified (dist_feedback_gain, concat_input_vector). The
  controller will be saturated to the sys.input_space using projection. 
 
  For an open-loop controller, only an concat_input_vector may be specified and
  dist_feedback_gain is set to zero. 
 
  See also examples/forwardStochasticReachCWH.m, examples/cwhSReachPointDemo.m
 
  =============================================================================
  [concat_state_realization, concat_disturb_realizations] = ...
     generateMonteCarloSims(n_monte_carlo_sims, sys, initial_state, ...
        time_horizon, optimal_input_vector, optimal_input_gain)
 
  Inputs:
  -------
    n_monte_carlo_sims  - Number of Monte-Carlo simulation particles to be used 
                          for estimation of the reach-avoid probability
    sys                 - System description as a LtiSystem/LtvSystem object
    initial_state       - Deterministic x_0
    time_horizon        - Time horizon (N) of the stochastic reach-avoid problem
    concat_input_vector - [Optional] Open-loop controller, a column vector of
                          dimension (sys.input_dim*N) x 1 | Required only if 
                          the system is controlled
    dist_feedback_gain  - [Optional] Affine disturbance feedback gain for the
                          concatenated disturbance vector, a matrix of dimension
                          (sys.input_dim*N) x (sys.dist_dim*N) | Required only
                          if the system is controlled, the controller is affine
                          disturbance feedback, and the gain matrix must be
                          lower block triangular (with zeros in its block
                          diagonal elements) for causality | See Notes
    verbose             - [Optional] Verbosity of this function when saturating
                          affine disturbance feedback controllers
 
  Outputs:
  --------
    concat_state_realization  - Matrix of concatenate state (row) vectors
                                stacked columnwise. Each row comprises of the
                                state trajectory as [x_0; x_1; x_2; ...; x_N]
    concat_disturb_realization- Matrix of concatenate disturbance (row) vectors
                                stacked columnwise. Each row comprises of
                                the state trajectory as [w_0; w_1; ...; w_{N-1}]
    saturation_indx           - Binary vector that indicates which realizations
                                had their associated affine disturbance feedback
                                controller saturated. Potentially non-zero only
                                if the input_gain is non-empty
 
  Notes:
  ------
  * Assumes IID disturbance for the LTI/LTV system. 
  * For controlled system, an open-loop controller NEEDS to be provided. The
    optimal_input_vector should be a ((sys.input_dim) * time_horizon)-dim.
    vector U = [u_0; u_1; ...; u_N] (column vector).
  * For uncontrolled system, the optimal_input_vector NEED NOT be provided
    dist_feedback_gain must be lower
  * The disturbance feedback gain matrix must be lower block triangular, with
    its block diagonal submatrices as zero matrices. This ensures that the
    affine disturbance feedback controller's value at any point of time depends
    only on the past disturbance values => causal controller. This function DOES
    NOT check for this structure in the input gain. TODO
  * Affine disturbance feedback controllers CAN NOT satisfy hard control bounds
    when the disturbance is unbounded (like Gaussian). Therefore, we will
    saturate the controller realization (associated with the disturbance
    realization) via projection on to the concatenated input space.
    Specifically, we solve the corresponding optimization problem for each
    concatenated disturbance realization W
    
        minimize || U - (MW + D)||_2
        subject to 
            U \in \mathcal{U}^T
 
    where U is the decision variable, \mathcal{U} is the input space, T is the
    time horizon, M is the affine disturbance feedback gain, and D is the affine
    disturbance feedback bias.
 
  ============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```
