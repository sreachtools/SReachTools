---
layout: docs
title: SReachSetCcO.m
---

```
  Compute the stochastic reach set corresponding to the stochastic reachability 
  problem of a target tube using convex chance-constraint optimization
  =============================================================================
 
  SReachSetCcO computes the open-loop controller-based underapproximative
  stochastic reach set to the problem of stochastic reachability of a target
  tube as discussed in
 
  A. Vinod and M. Oishi, "Scalable underapproximative verification of stochastic
  LTI systems using convexity and compactness," In Proc. Hybrid Syst.: Comput. &
  Ctrl., pages 1--10, 2018. HSCC 2018
 
  A. Vinod and M. Oishi, "Stochastic reachability of a target tube: Theory and
  computation," IEEE Transactions in Automatic Control, 2018 (submitted)
  https://arxiv.org/pdf/1810.05217.pdf.
 
  =============================================================================
 
  [polytope, extra_info] = SReachSetCcO(method_str, sys, prob_thresh,...
    safety_tube, options)
  
  Inputs:
  -------
    method_str  - Solution technique to be used. Must be 'chance-open'
    sys         - System description (LtvSystem/LtiSystem object)
    prob_thresh - Probability threshold at which the set is to be constructed
    safety_tube - Collection of (potentially time-varying) safe sets that
                  define the safe states (Tube object)
    options     - Collection of user-specified options for 'chance-open'
                  (Matlab struct created using SReachSetOptions)
 
  Outputs:
  --------
    polytope   - Underapproximative polytope of dimension sys.state_dim which
                 underapproximates the stochastic reach set
    extra_info - A list of Matlab structs that comprises of auxillary
                 information from the set computation.
                 The list has two members extra_info_wmax and extra_info_cheby.
                 Both of these structs contain the following information:
                    1. xmax - Initial state that has the maximum reach
                              probability to stay with the safety tube using an
                              open-loop controller (via the method in use)
                    2. Umax - Optimal open-loop policy ((sys.input_dim) *
                              time_horizon)-dimensional vector 
                              U = [u_0; u_1;...; u_N] (column vector) for xmax
                              (via the method in use)
                    3. xmax_reach_prob 
                            - Maximum attainable reach probability to
                              stay with the safety tube using an open-loop
                              controller
                    4. opt_theta_i 
                            - Vector comprising of scaling factors along each
                              user-specified direction of interest
                    5. opt_input_vec_at_vertices 
                            - Optimal open-loop policy ((sys.input_dim) *
                              time_horizon)-dim.  vector U = [u_0; u_1; ...;
                              u_N] (column vector) for each vertex of the
                              polytope
                    6. opt_reach_prob_i
                            - Maximum attainable reach probability to stay with
                              the safety tube at the vertices of the polytope
                    7. vertices_underapprox_polytope
                            - Vertices of the polytope
                                xmax + opt_theta_i * options.set_of_dir_vecs
 
  Notes:
  ------
  * extra_info(1).xmax_reach_prob is the highest prob_thresh that may be given
    while obtaining a non-trivial underapproximation
  * See @LtiSystem/getConcatMats for more information about the
      notation used.
  
  =============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
 
```
