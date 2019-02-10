---
layout: docs
title: spreadPointsOnUnitSphere.m
---

```
  Computes a collection of n_points vectors of dimension n_dim that have a large
  minimum pairwise-separation
  ============================================================================
  
  Given n_points, we spread n_points_first_quad = ceil(n_points - 2*n_dim) in
  the first quadrant using the following non-convex optimization problem,
  
  maximize r
  subject to
                   r >= 0,                                                   (1)
   || x_i - x_j ||_2 >= r,   i,j\in{1,..,n_points_first_quad}, i<j           (2) 
   || x_i - e_j ||_2 >= r,   i\in{1,..,n_points_first_quad}, j\in{1,...n_dim}(3)   
         || x_i ||_2 <= 1,   i\in{1,..,n_points_first_quad}                  (4)  
         || x_i ||_2 >= 0.8, i\in{1,..,n_points_first_quad}                  (5)   
                 x_i >= r/2                                                  (6)  
  where e_j refers to the standard vector (zeros with 1 at jth position). Here,
  (1) enforces positive separation, (2) and (3) enforces the smallest
  pairwise separation is above r (among each other and the standard vectors),
  (4) and (5) approximates || x_i ||_2 = 1, and (6) enforces the separation
  constraint is satisfied even among the reflections/rotations.
 
  This optimization problem is non-convex, and we solve it to a local optimality
  using difference-of-convex approach. Specifically, constraints (2), (3), and
  (5), which are reverse-convex, are tightened to their first-order Taylor
  series (under)approximation and the resulting linear constraints are enforced
  in their place. This method is discussed in:
 
       J. D. Gleason, A. P. Vinod, and M. M. K. Oishi. 2018. Lagrangian 
       Approximations for Stochastic Reachability of a Target Tube. 
       online. (2018). https://arxiv.org/abs/1810.07118
 
  Next, we reflect/rotate these vectors in the first quadrant to occupy in all
  other quadrants. In the end, we tack on e_j and -e_j for each dimension j
  (hence, the - 2*n_dim).
 
  ============================================================================
  
  [opt_locations, separation] = spreadPointsOnUnitSphere(n_dim,n_points,verbose)
 
  Inputs:
  -------
  n_dim     - Dimension of the unit sphere on which we wish to spread the points 
  n_points  - Number of points we wish to spread the points (Will be rounded up
              to the smallest k such that 2^(n_dim) k + 2*n_dim >= n_points
  verbose   - Verbosity of this function
 
  Outputs:
  --------
  opt_locations - Unit vectors given as a n_dim x n_points
  separation    - The minimum pairwise-separation across the vectors
 
  Notes:
  ------
  * We enforce x_i^T x_i >= 0.8^2 instead of x_i^T x_i >= 1, so that the problem
    converges faster.
  ============================================================================
 
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
 
```
