---
layout: docs
title: computeNormCdfInvOverApprox.m
---

```
  Compute a piecewise-linear overapproximation of norminv(1-x) for 
  x \in [lb_delta,0.5] to a user-specified quality (Internal function)
  =============================================================================
 
  computeNormCdfInvOverApprox generates a piecewise-linear overapproximation of
  norminv(1-x) for x\in[lb_delta,0.5]. Specifically, given any
  z\in[lb_delta,0.5],
 
  norminv(1-x) + pwa_accuracy>max(cdf_approx_m * x + cdf_approx_c)>norminv(1-x),
 
  with lb_delta = max(max_delta/n_lin_consts/10,1e-8).
 
  This function implements Algorithm 2 of the following paper:
 
  A. Vinod and M. Oishi. Affine controller synthesis for stochastic reachability
  via difference of convex programming. In Proc. Conf. Dec. & Ctrl. 2019.
  (submitted). https://hscl.unm.edu/affinecontrollersynthesis/
 
  =============================================================================
 
  
  Inputs: 
  -------
    max_delta    - This is the maximum tolerance for violation of the joint
                   chance constraint | risk allocation can not exceed this value
    pwa_accuracy - Accuracy expected from the chance constraint formulation with 
                   n_lin_consts, no. of individual chance constraints
    n_lin_consts - No. of individual chance constraints
 
  Outputs:
  --------
    overapprox_m - Secant slopes that will overapproximate norminv(1-x)
    overapprox_c - Secant y-intercepts that will overapproximate norminv(1-x)
    lb_phiinv    - Lower bound on x in norminv(1-x) for which the provided PWA
                    approximation
    norminv_knots- Breakpoints of the PWA overapproximation, i.e., points
                    at which the PWA overapproximation coincides with the 
                    norminv curve
 
  Notes:
  * Partial INPUT HANDLING: Checks only if the bounds provided are accurate
  * MATLAB DEPENDENCY: None
  * The requested desired_accuracy/n_lin_consts/10 can not be smaller than 1e-8.
  
  ==============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
 
```
