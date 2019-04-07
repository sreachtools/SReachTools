---
layout: docs
title: SReachPointVoO.m
---

```
  Solve the problem of stochastic reachability of a target tube (a lower bound
  on the maximal reach probability and an open-loop controller synthesis) using
  undersampled particle filter control
  =============================================================================
 
  SReachPointVoO implements a mixed-integer linear program-based approximation
  to the stochastic reachability of a target tube problem. To improve the
  tractability, we undersample the particles by utilizing the kmeans
  clustering/Voronoi partitioning of the disturbance samples. 
 
  H. Sartipizadeh, A. Vinod,  B. Acikmese, and M. Oishi, "Voronoi
  Partition-based Scenario Reduction for Fast Sampling-based Stochastic
  Reachability Computation of LTI Systems", In Proc. Amer. Ctrl. Conf., 2019
  
  In contrast to `particle-open' approach, implemented in SReachPointPaO and
  described in
 
  K. Lesser, M. Oishi, and R. Erwin, "Stochastic reachability for control of
  spacecraft relative motion," in IEEE Conference on Decision and Control (CDC),
  2013,
 
  SReachPointVoO computes an underapproximation of the maximal reach probability
  with the probabilistic guarantee that the risk of failure of this
  underapproximation is no larger than the user-specified failure_risk. 
 
     High-level desc.   : Sample scenarios based on the additive noise and solve
                          a mixed-integer linear program to make the maximum
                          number of scenarios satisfy the reachability
                          objective.  In addition, we use Voronoi partition to
                          drastically improve the tractability while preserving
                          the underapproximation quality
     Approximation      : Overapproximation bounded above by a user-specified
                          tolerance
     Controller type    : Open-loop controller that satisfies the hard input
                          bounds
     Optimality         : Optimal (w.r.t scenarios drawn) open-loop controller
                          for the underapproximation problem 
 
  =============================================================================
 
  [approx_stoch_reach, opt_input_vec, varargout] = SReachPointVoO(sys, ...
    initial_state, safety_tube, options)
 
  Inputs:
  -------
    sys          - System description (LtvSystem/LtiSystem object)
    initial_state- Initial state for which the maximal reach probability must be
                   evaluated (A numeric vector of dimension sys.state_dim)
    safety_tube  - Collection of (potentially time-varying) safe sets that
                   define the safe states (Tube object)
    options      - Collection of user-specified options for 'voronoi-open'
                   (Matlab struct created using SReachPointOptions)
 
  Outputs:
  --------
    approx_stoch_reach 
                - An approximation of the stochastic reachability of a target
                  tube problem computed using undersampled particle control
                  approach using kmeans. In contrast to `particle-open'
                  approach, this approximation permits a user-defined
                  upper-bound on the overapproximation error. While it is
                  expected to lie in [0,1], it is set to -1 in cases where the
                  CVX optimization fails (cvx_status \not\in {Solved,
                  Inaccurate/Solved}) or if the optimal value is below
                  max_overapprox_err.  
    opt_input_vec
                - Open-loop controller: column vector of dimension
                  (sys.input_dim*N) x 1
    extra_info  - A MATLAB struct containing the information about partitioning
                  of GW space. The struct contains the following info:
                   n_particles - Number of particles based off Hoeffding's
                                 inequality
                   n_kmeans    - Number of bins for kmeans clustering
                   GW_centroids- Centroids obtained from kmeans clustering
                   GW_realizations
                               - Realizations for the random vector GW
                   kmeans_idx  - Output of MATLAB's kmeans function | Index
                                 of the Voronoi centroids for each of the
                                 samples
                   approx_voronoi_stoch_reach
                               - Optimal solution of the undersampled MILP,
                                 before the reevaluation of the policy for
                                 tightening
                   offline_compute_time
                               - Time taken for computing the minimum required 
                                 samples, generating the samples, and computing
                                 the Voronoi centers and buffers (the steps
                                 that can be done independent of the
                                 initial state or the open-loop control
                                 vector => offline computable steps)
                   online_compute_time
                               - Time taken for solving the MILP and computing
                                 the tightened underapproximation
 
  See also SReachPoint.
 
  Notes:
  * This function requires CVX with Gurobi as the backend solver for optimizing
    the resulting mixed-integer linear program.
  * This function requires kmeans function which is part of MATLAB's
    Statistical and Machine Learning toolbox.
  * See @LtiSystem/getConcatMats for more information about the notation used.
  * This function computes an approximate stochastic reach probability with a
    probabilistic guarantee of being an underapproximation. Specifically, it
    computes an underapproximation of the maximal reach probability with the
    probabilistic guarantee that the risk of failure of this underapproximation
    is no larger than the user-specified failure_risk. 
  * The number of scenarios required to guarantee the above statement is:
        - directly proportional to the -log(failure_risk). 
        - inversely proportional to the square of maximum overapproximation
          error.
    The maximum overapproximation error is the correction term, the threshold
    that is removed from the optimal solution of the MILP to guarantee the
    underapproximation.
  * Ideally, we should solve the MILP with the number of scenarios prescribed by
    the above note. However, this number is typically high rendering the MILP
    intractable.
  * Therefore, we solve an undersampled MILP, where the user specifies the
    number of representative scenarios to consider (n_kmeans). Larger n_kmeans
    implies higher accuracy for the approx_stoch_reach. However, this comes at
    the computational cost of solving the MILP, which has n_kmeans binary
    variables.
  
  ============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```
