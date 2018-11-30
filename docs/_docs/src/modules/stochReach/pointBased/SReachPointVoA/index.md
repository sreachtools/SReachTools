---
layout: docs
title: SReachPointVoA.m
---

```
  Solve the problem of stochastic reachability of a target tube (a lower bound
  on the maximal reach probability and an affine controller synthesis) using
  under-sampled particle control and coordinate descent algorithm
  =============================================================================
 
  SReachPointVoA implements the under-sampled particle control to the problem of
  stochastic reachability of a target tube to construct an affine controller.
 
     High-level desc.   : Sample particles based on the additive noise and solve
                          a mixed-integer linear program to make the maximum
                          number of particles satisfy the reachability objective
                          In addition, we use Voronoi partition to
                          drastically improve the tractability while
                          preserving the underapproximation quality
     Approximation      : Overapproximation bounded above (in probability) by a
                          user-specified tolerance
     Controller type    : A history-dependent affine controller that satisfies
                          softened input constraints (controller satisfies the
                          hard input bounds upto a user-specified probabilistic
                          threshold)
     Optimality         : Suboptimal (w.r.t particles drawn) affine disturbance
                          feedback controller 
     Dependency (EXT)   : CVX, Gurobi
     SReachTool function: SReachPointVoA
     Paper              : TODO
 
  =============================================================================
 
    [lb_stoch_reach, opt_input_vec, opt_input_gain, ...
     risk_alloc_state, risk_alloc_input] = SReachPointVoA(sys,...
        initial_state, safety_tube, options)
  
  Inputs:
  -------
    sys          - System description (LtvSystem/LtiSystem object)
    initial_state- Initial state for which the maximal reach probability must be
                   evaluated (A numeric vector of dimension sys.state_dim)
    safety_tube  - Collection of (potentially time-varying) safe sets that
                   define the safe states (Tube object)
    options      - Collection of user-specified options for 'voronoi-affine'
                   (Matlab struct created using SReachPointOptions)
 
  Outputs:
  --------
    lb_stoch_reach 
                - Lower bound on the stochastic reachability of a target tube
                  problem computed using chance constraints and
                  difference-of-convex techniques
    opt_input_vec, 
      opt_input_gain
                - Controller U=MW+d for a concatenated input vector 
                    U = [u_0; u_1; ...; u_{N-1}] and concatenated disturbance
                    vector W=[w_0; w_1; ...; w_{N-1}]. 
                    - opt_input_gain: Affine controller gain matrix of dimension
                        (sys.input_dim*N) x (sys.dist_dim*N)
                    - opt_input_vec: Open-loop controller: column vector 
                      dimension
                        (sys.input_dim*N) x 1
    kmeans_info - A MATLAB struct containing the information about partitioning
                  of W space. The struct contains the following info:
                   n_particles    - Number of particles based off Hoeffding's
                                    inequality
                   n_kmeans       - Number of bins for kmeans clustering
                   W_centroids    - Centroids obtained from kmeans clustering
                   W_realizations - Realizations for the random vector W
 
 
  See also SReachPoint.
 
  Notes:
  * We recommend using this function through SReachPoint.
  * This function requires CVX to work.
  * This function returns a **lower bound to the maximal reach probability under
    hard input constraints**. This lower bound is obtained by a linear
    transformation of the maximal reach probability associated with the
    unsaturated affine controller using the user-specified likelihood threshold
    on the hard input constraints. See Theorem 1 of the paper cited above.
  * Due to numerical issues, we add a small positive perturbation to the
    b term, whenever determining containment --- Ax<=b
  * See @LtiSystem/getConcatMats for more information about the notation used.
  
  ============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```
