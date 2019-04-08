---
layout: docs
title: SReachPointOptions.m
---

```
  Create user-specifiable options for use with SReachPoint()
  =============================================================================
 
  SReachPointOptions creates a MATLAB struct that contains user-specifiable
  options that may be used with SReachPoint
 
  =============================================================================
 
    options = SReachPointOptions(prob_str, method_str, varargin)
  
  Inputs:
  -------
    prob_str    - String specifying the problem of interest. For each case, we
                  compute the optimal value function that maps initial states
                  to different maximal reach probabilities
                      1. 'term' : Stay within the safety_tube
    method_str  - Solution technique to be used (user-specifiable
                  options associated with each technique is enumerated)
                      'chance-open'  -- Convex chance-constrained approach for
                                        an open-loop controller synthesis
                                        1. pwa_accuracy: Accuracy of the
                                                piecewise affine approximation
                                                of norminvcdf used [Default:
                                                1e-3]
                      'genzps-open'  -- Genz's algorithm + Patternsearch
                                        1. desired_accuracy: Accuracy of
                                                Gaussian integral => Accuracy of
                                                the result [Default: 5e-2]
                                        2. PSoptions: MATLAB struct generated
                                                using psoptimset()
                                        3. thresh: An upper bound on useful
                                                reach probability. This can be 
                                                used to specify reach_prob >= 
                                                thresh \in (0,1] See notes 
                                                [Default: 1]
                      'particle-open'-- Particle control approach that uses
                                        mixed-integer linear programming
                                        1. n_particles: Number of particles to
                                                use [Default: 100] | Must
                                                be less than max_particles.
                                        2. bigM: A large positive constant value
                                                that is used in the mixed
                                                integer formulation [Default:
                                                5000]
                                        3. verbose: Verbosity of the 
                                                implementation (feedback for the
                                                user) | Takes values from 0 to 2
                                                [Default: 0]
                                        4. max_particles: Maximum particles
                                                permitted. This bound is
                                                used to throw a pre-emptive
                                                error for very demanding
                                                problem requirements
                                                [Default: 200]
                      'voronoi-open' -- Voronoi-based undersampling of particle
                                        control approach to compute open loop
                                        1. failure_risk: Risk of the
                                                probabilistic overapproximation
                                                bound failing [Default: 1e-4]
                                        2. max_overapprox_err: Maximum
                                                overapproximation error
                                                (probabilistically) tolerable up
                                                to the failure_risk
                                                [Default: 1e-2]
                                        3. n_kmeans: Number of kmeans cluster
                                                points/ Voronoi centers
                                                [Default: 30]
                                        4. bigM: A large positive constant value
                                                used in the mixed integer 
                                                formulation [Default: 100]
                                        5. verbose: Verbosity of the 
                                                implementation (feedback for the
                                                user) | Takes values from 0 to 2
                                                [Default: 0]
                                        6. max_particles: Maximum particles
                                                permitted. This bound is
                                                used to throw a pre-emptive
                                                error for very demanding
                                                problem requirements
                                                [Default: 1e5]
                      'chance-affine'-- Difference-of-convex chance-constrained 
                                        approach for affine controller synthesis
                                        1. [MUST HAVE] max_input_viol_prob:
                                                Probabilistic relaxation of the
                                                hard input constraints 
                                                [Default: 1e-2]
                                        2. verbose: Verbosity of the 
                                                implementation (feedback for the
                                                user) | Takes values from 0 to 2
                                                [Default: 0]
                                        3. pwa_accuracy: Accuracy of the
                                                piecewise affine approximation
                                                of norminvcdf used [Default:
                                                1e-3]
                                        Difference-of-convex parameters: 
                                        4. tau_initial: Initialization of the 
                                                slack multiplier [Default: 1]
                                        5. scaling_tau: Scaling factor to the 
                                                slack multiplier [Default: 2]
                                        6. tau_max: Maximum value for the 
                                                scaling factor [Default: 1e5]
                                        7. iter_max: Maximum number of
                                                iterations for the difference of
                                                convex iterative algorithm 
                                                [Default: 200]
                                        8. dc_conv_tol: Tolerance for exiting 
                                                the iterative algorithm
                                                [Default: 1e-4]
                                        9. slack_tol: Tolerance for the sum
                                                of slack vars for penalty DC
                                                [Default: 1e-8]
                      'chance-affine-uni'
                                     -- Uniform risk allocation for affine 
                                        controller synthesis
                                        1. [MUST HAVE] max_input_viol_prob:
                                                Probabilistic relaxation of the
                                                hard input constraints 
                                                [Default: 1e-2]
                                        2. verbose: Verbosity of the 
                                                implementation (feedback for the
                                                user) | Takes values from 0 to 1
                                                [Default: 0]
                                        3. state_bisect_tol: Bisection
                                                tolerance for the state 
                                                constraint [Default: 1e-3]
                                        4. input_bisect_tol: Bisection
                                                tolerance for the input 
                                                constraint [Default: 1e-3]
 
  Outputs:
  --------
    options     - Collection of user-specified options for given method_str
 
  See also SReachPoint.
 
  Notes:
  * SReachPoint() will call this function internally using the default
    values if SReachPointOptions()-based options is not explicitly provided
    to SReachPoint().
  * To specify a desired set of samples V to use when undersampling in 
    voronoi-X, set the undersampling fraction to be very small (say 1e-4/1e-5) 
    and set min_samples to V.
  * For sampling-based approaches, we impose a heuristic maximum allowed
    particles. This bound is used to throw a pre-emptive error for very 
    demanding problem requirements. The user MAY MODIFY it to allow the
    computations beyond the limit.
  * max_input_viol_prob has been left as addParameter instead of
    addRequired, so that default values may be specified.
  * For 'genzps-open', options.thresh permits early termination if instead
    of maximal reach probability is of interest. While this function is for 
    maximal reach probability, this optional feature permits the reuse of
    code in SReachSet.
  ============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```
