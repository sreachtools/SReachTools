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
            <li class="doc-list"><a href="#SReachEllipsoid-method-mtimes">mtimes</a></li>
            <li class="doc-list"><a href="#SReachEllipsoid-method-support_fun">support_fun</a></li>
        </ul>
    </ul>
</ul>

{:#SReachEllipsoid}
### SReachEllipsoid
```
  Creates an ellipsoid (x - c)^T Q^{-1} (x-c) <= 1
 
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

### Method: support_fun
{:#SReachEllipsoid-method-support_fun}
```
  Support function of the ellipsoid object
  ====================================================================
 
  Inputs:
  -------
    l  - A query column vector or a collection of query vectors stacked 
         as rows
 
  Outputs:
  --------
    support_fun_val - max_{y \in ellipsoid} l'*y
 
  =====================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
  
```

