---
layout: docs
title: SReachEllipsoid.m
---

<ul class="doc-list">
    <li class="doc-list"><a href="#SReachEllipsoid">SReachEllipsoid</a></li>
    <ul class="doc-list">
        <li><a href="#SReachEllipsoid-SReachEllipsoid">Constructor</a></li>
        <li>Properties</li>
        <ul class="doc-list">
            <li class="doc-list"><a href="#SReachEllipsoid-prop-center">center</a></li>
            <li class="doc-list"><a href="#SReachEllipsoid-prop-shape_matrix">shape_matrix</a></li>
            <li class="doc-list"><a href="#SReachEllipsoid-prop-dim">dim</a></li>
        </ul>
        <li>Methods</li>
        <ul class="doc-list">
            <li class="doc-list"><a href="#SReachEllipsoid-method-contains">contains</a></li>
            <li class="doc-list"><a href="#SReachEllipsoid-method-plus">plus</a></li>
            <li class="doc-list"><a href="#SReachEllipsoid-method-mtimes">mtimes</a></li>
            <li class="doc-list"><a href="#SReachEllipsoid-method-support">support</a></li>
        </ul>
    </ul>
</ul>

{:#SReachEllipsoid}
### SReachEllipsoid
```
  Creates an ellipsoid (x - c)^T Q^{-1} (x-c) <= 1
  =============================================================================
  
  Create an ellipsoid object defined by the equation
  
    E = { x \in R^{n} : (x - c)^{T} Q^{-1} (x - c) <= 1 }
  
  These ellipsoid objects are often used when creating bounded disturbances
  for the SReachSet Lagrangian methods when using Gaussian disturbances
  
  Usage:
  ------
  % Create unit ellipsoid
  sre = SReachEllipsoid([0;0], eye(2));
  
  =============================================================================
 
  SReachEllipsoid Properties:
  ---------------------------
    center                          - Center of the ellipsoid (c)
    shape_matrix                    - Shape matrix of the ellipsoid Q
    dim                             - Dimension of the ellipsoid
  
  SReachEllipsoid Methods:
  ------------------------
    SReachEllipsoid/SReachEllipsoid - Constructor
    support                         - Support function of the ellipsoid
    contains                        - Checks if a point (column vector) or a
                                      collection of points (matrix of column
                                      vectors) is within the ellipsoid
 
  Apart from these methods, the following commands work
    disp                            - Displays critical info about the ellipsoid
    F * ell | ell * F               - Multiplication of ellipsoid by a n x dim -
                                      dimensional matrix or a scalar F
    F + ell | ell + F | ell + poly  - Add a deterministic vector/scalar to an
                                      ellipsoid | Overapproximate the Minkowski
                                      sum of a polyhedron with ellipsoid
  
  Notes:
  ------
  * The ellipsoid can be full-dimensional (Q non-singular) or be a 
    lower-dimensional ellipsoid embedded in a high dimensional space (Q 
    singular)
 
  ===========================================================================
 
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
  

    Reference page in Doc Center
       doc SReachEllipsoid

```

{:#SReachEllipsoid-SReachEllipsoid}
### Constructor
```
  Creates an ellipsoid (x - c)^T Q^{-1} (x-c) <= 1
  =============================================================================
  
  Create an ellipsoid object defined by the equation
  
    E = { x \in R^{n} : (x - c)^{T} Q^{-1} (x - c) <= 1 }
  
  These ellipsoid objects are often used when creating bounded disturbances
  for the SReachSet Lagrangian methods when using Gaussian disturbances
  
  Usage:
  ------
  % Create unit ellipsoid
  sre = SReachEllipsoid([0;0], eye(2));
  
  =============================================================================
 
  SReachEllipsoid Properties:
  ---------------------------
    center                          - Center of the ellipsoid (c)
    shape_matrix                    - Shape matrix of the ellipsoid Q
    dim                             - Dimension of the ellipsoid
  
  SReachEllipsoid Methods:
  ------------------------
    SReachEllipsoid/SReachEllipsoid - Constructor
    support                         - Support function of the ellipsoid
    contains                        - Checks if a point (column vector) or a
                                      collection of points (matrix of column
                                      vectors) is within the ellipsoid
 
  Apart from these methods, the following commands work
    disp                            - Displays critical info about the ellipsoid
    F * ell | ell * F               - Multiplication of ellipsoid by a n x dim -
                                      dimensional matrix or a scalar F
    F + ell | ell + F | ell + poly  - Add a deterministic vector/scalar to an
                                      ellipsoid | Overapproximate the Minkowski
                                      sum of a polyhedron with ellipsoid
  
  Notes:
  ------
  * The ellipsoid can be full-dimensional (Q non-singular) or be a 
    lower-dimensional ellipsoid embedded in a high dimensional space (Q 
    singular)
 
  ===========================================================================
 
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
  

    Reference page in Doc Center
       doc SReachEllipsoid

```

### Property: center
{:#SReachEllipsoid-prop-center}
```
  SReachEllipsoid/center
  ==================================================================
  
  Column vector indicating the center of the ellipsoid
 
```

### Property: shape_matrix
{:#SReachEllipsoid-prop-shape_matrix}
```
  SReachEllipsoid/shape_matrix
  ==================================================================
  
  Shape matrix for the ellipsoid
  
```

### Property: dim
{:#SReachEllipsoid-prop-dim}
```
  SReachEllipsoid/dim
  ==================================================================
  
  Dimension of the ellipsoid dimension
  
```

### Method: contains
{:#SReachEllipsoid-method-contains}
```
  Checks if a point (column vector) or a collection of points (matrix of
  column vectors) is within the ellipsoid
  ====================================================================
  
  Inputs:
  -------
    obj         - Ellipsoid object
    test_points - Point (column vector) or a N collection of points
                  (matrix of column vectors) is within the ellipsoid
 
  Outputs:
  --------
    newobj      - Boolean vector Nx1 that describe the containment
 
  Notes:
  ------
  * Requires CVX for vectorized norm.
 
  ====================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```

### Method: plus
{:#SReachEllipsoid-method-plus}
```
  Override of MATLAB plus command
  ====================================================================
  
  Inputs:
  -------
    obj - Ellipsoid object
    v   - Deterministic vector to be added to the random vector OR
          a Polytope object
 
  Outputs:
  --------
    newobj - Ellipsoid obj (obj + v) for deterministic vector/scalar v
             Polyhedron obj (obj \oplus v) for polytopic v (overapprox)
 
  Notes:
  ------
  * For a polytopic v, newobj is an (Polyhedron overapproximation of the
    minkowski sum, computed via sampling the support function.
 
  ====================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```

### Method: mtimes
{:#SReachEllipsoid-method-mtimes}
```
  Override of MATLAB multiplication command
  ====================================================================
  
  Inputs:
  -------
    obj - SReachEllipsoid object
    F   - Linear transformation matrix for multiplication
 
  Outputs:
  --------
    newobj - SReachEllipsoid object (F*obj)
 
  ====================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```

### Method: support
{:#SReachEllipsoid-method-support}
```
  Support function of the ellipsoid object
  ====================================================================
 
  Inputs:
  -------
    l   - A query column vector or a collection of query vectors stacked 
          as columns
 
  Outputs:
  --------
    val - max_{y \in ellipsoid} l'*y
 
  =====================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
  
```

