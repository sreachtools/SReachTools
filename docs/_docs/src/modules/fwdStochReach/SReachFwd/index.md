---
layout: docs
title: SReachFwd.m
---

```
  Perform forward stochastic reachability analysis of a linear system with
  additive stochastic disturbance, initialized to a random vector
  ============================================================================
 
  Forward stochastic reachability analysis quantifies the stochasticity of the
  state at a future time instant, as well as the probability of the state lying
  in a pre-specified set at a future time instant. This function implements
  ideas from
 
  A. Vinod, B. HomChaudhuri, and M. Oishi, "Forward Stochastic Reachability
  Analysis for Uncontrolled Linear Systems using Fourier Transforms", In
  Proceedings of the 20th International Conference on Hybrid Systems:
  Computation and Control (HSCC), 2017.
 
  While the Fourier transform-based results apply for arbitrary 
  distributions, SReachFwd considers only Gaussian-perturbed LTI systems. 
  In this case, the approach coincides with Kalman filter updates, and it 
  is grid-free, recursion-free, and sampling-free. 
 
  SReachFwd also exploits the functionality of random vector to provide
  forward stochastic reachability analysis of linear systems perturbed by
  non-Gaussian disturbance. In this case, we use Monte-Carlo simulations to
  estimate the mean, covariance, and the probability of lying in a set.
 
  See also examples/cwhSReachFwd.m.
 
  ============================================================================
  
  varargout = SReachFwd(prob_str, sys, initial_state, target_time, varargin)
  
  Inputs:
  -------
    prob_str      - String specifying the problem of interest
                        1. 'state-stoch' : Provide mean and covariance at state
                                           at the spjecified time
                        2. 'state-prob'  : Compute the probability that the
                                           state will lie in a polytope at the
                                           specified time
                        3. 'concat-stoch': Provide mean and covariance of the
                                           concatenated state vector up to a
                                           specified time
                        4. 'concat-prob' : Compute the probability that the
                                           concatenated state vector up to a
                                           specified time lies in the given
                                           target tube
    sys           - System description as a LtiSystem/LtvSystem object
    initial_state - Initial state as a deterministic n-dimensional vector
                    or a RandomVector object
    target_time   - Time of interest (non-negative scalar) | If target_time = 0,
                    then the stochasticity of the initial state is analyzed.
    target_set/tube  
                  - [Required only for state/concat-prob] Polyhedron/Tube object
                    over which the probability must be computed
    desired_accuracy 
                  - [Optional for state/concat-prob] Maximum absolute deviation
                    from the true probability estimate [Default: 1e-2]
 
  Outputs:
  --------
    rv            - ['state/concat-stoch'] Random vector describing the
                    state / concatenated state vector | It is a vector of
                    dimension (target_time + 1) * sys.state_dim to include the
                    initial state
    prob          - ['state/concat-prob'] Probability of occurence
 
  Notes:
  ------
  * Assumes IID disturbance.
  * The outputs are either (mean_vec, cov_mat) or (prob), depending on the
    method_str
  * For concat-prob, target_time can be any value between 0 and N where
    target_tube has a length of N+1 (N+1 target sets to include the
    constraints on the initial state).
  * For XXX-prob, the random vector is provided as the second output.
  * Prob{ w \in \theta Polytope(A,b) } is computed using
    RandomVector/getProbPolyhedron.
  * Requires the desired_accuracy to be at least 1e-2.
 
  ============================================================================
 
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
 
```
