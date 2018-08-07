---
layout: docs
title: computeNormCdfInvOverApprox.m
---

```
  SReachTools/stochasticReachAvoid/computeNormCdfInvOverApprox: Compute a
  piecewise-linear overapproximation of norminv(1-x) for x\in[1e-5,0.5] to the
  quality of 1e-4
  =============================================================================
 
  computeNormCdfInvOverApprox generates a piecewise-linear overapproximation of
  norminv(1-x) for x\in[1e-5,0.5]. Specifically, given any z\in[1e-5,0.5],
  norm(1-z) +1e-4 > max(cdf_approx_m * z + cdf_approx_c) > norm(1-z).
 
  USAGE: See getUnderapproxStochReachAvoidSet,
  computeCcLowerBoundStochReachAvoidPwlRisk.
 
  =============================================================================
 
  
  Inputs: None
  -------
 
  Outputs:
  --------
    cdf_approx_m - Secant slopes that will overapproximate norminv(1-x)
    cdf_approx_c - Secant y-intercepts that will overapproximate norminv(1-x)
    lb_x         - (Optional) Lower-bound on the range of the approximation
    diff_val     - (Optional) Smallest step-size between the secant end points
 
  Notes:
  * NOT ACTIVELY TESTED: TODO
  * NO INPUT HANDLING: No arguments needed
  * MATLAB DEPENDENCY: None
  * Max error of approximation is 0.99818e-04 (estimated to a nbd of 1e-7)
  * The end points of the secants are obtained by the sequence 
        {lb_x, lb_x + h gamma^{0:n_x},0.5, and the midpoints} 
    where lb_x =1e-5, h=1e-6, gamma=1.088, n_x =|_log((0.5-lb_x)/h)/log(gamma)_|
  * If we desire to optimize delta_i where i\in [1,M], then this approach
    results in a maximum artificial conservativeness (on top of Boolean
    conservativeness) of 1e-4*M.
  
  =============================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
 
```
