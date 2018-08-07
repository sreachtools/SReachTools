---
layout: docs
title: SimpleBox.m
---

<ul class="doc-list">
    <li class="doc-list"><a href="#SimpleBox">SimpleBox</a></li>
    <ul class="doc-list">
        <li><a href="#SimpleBox-SimpleBox">Constructor</a></li>
        <li>Properties</li>
        <ul class="doc-list">
            <li class="doc-list"><a href="#SimpleBox-prop-vertices">vertices</a></li>
            <li class="doc-list"><a href="#SimpleBox-prop-center">center</a></li>
            <li class="doc-list"><a href="#SimpleBox-prop-dx">dx</a></li>
            <li class="doc-list"><a href="#SimpleBox-prop-dim">dim</a></li>
        </ul>
        <li>Methods</li>
        <ul class="doc-list">
            <li class="doc-list"><a href="#SimpleBox-method-getPolyhedron">getPolyhedron</a></li>
            <li class="doc-list"><a href="#SimpleBox-method-computeGaussianProbability">computeGaussianProbability</a></li>
            <li class="doc-list"><a href="#SimpleBox-method-getBounds">getBounds</a></li>
        </ul>
    </ul>
</ul>

{:#SimpleBox}
### SimpleBox
```
  SReachTools/SimpleBox  Class definition to obtain vertices of a n-dimensional 
  wbox
  ===========================================================================
 
  Class to obtain vertices of an n-dimensional box; often used for computing
  probabilities in dynamic programming recursions
 
  Usage
  -----
  % call by passing in vertices
  simpbox = SIMPLEBOX([1,1;-1,-1;-1,1;1,-1]);
  
  % call by passing center and deltas 
  simpbox = SIMPLEBOX(0, [1, 1])
  
  ===========================================================================
 
  SIMPLEBOX Properties:
  ---------------------
    vertices- Array (m x n) of vertices; each vertex is a (1 x n) array
    center  - Array (1 x n) of box center location
    dx      - Array (1 x n) of half-lengths of box sides
    dim     - Dimension of box (scalar)
 
  SIMPLEBOX Methods:
  ------------------
    SimpleBox/SimpleBox        - Class constructor
    getPolyhedron              - Get Polyhedron object for box
    computeGaussianProbability - Compute the probability of Gaussian random 
                                 variable being in box
 
  ===========================================================================
 
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
  

    Reference page in Doc Center
       doc SimpleBox

```

{:#SimpleBox-SimpleBox}
### Constructor
```
  SReachTools/SimpleBox  Class definition to obtain vertices of a n-dimensional 
  wbox
  ===========================================================================
 
  Class to obtain vertices of an n-dimensional box; often used for computing
  probabilities in dynamic programming recursions
 
  Usage
  -----
  % call by passing in vertices
  simpbox = SIMPLEBOX([1,1;-1,-1;-1,1;1,-1]);
  
  % call by passing center and deltas 
  simpbox = SIMPLEBOX(0, [1, 1])
  
  ===========================================================================
 
  SIMPLEBOX Properties:
  ---------------------
    vertices- Array (m x n) of vertices; each vertex is a (1 x n) array
    center  - Array (1 x n) of box center location
    dx      - Array (1 x n) of half-lengths of box sides
    dim     - Dimension of box (scalar)
 
  SIMPLEBOX Methods:
  ------------------
    SimpleBox/SimpleBox        - Class constructor
    getPolyhedron              - Get Polyhedron object for box
    computeGaussianProbability - Compute the probability of Gaussian random 
                                 variable being in box
 
  ===========================================================================
 
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
  

    Reference page in Doc Center
       doc SimpleBox

```

### Property: vertices
{:#SimpleBox-prop-vertices}
{% include important-note.html content="Property currently has no help documentation." %}

### Property: center
{:#SimpleBox-prop-center}
{% include important-note.html content="Property currently has no help documentation." %}

### Property: dx
{:#SimpleBox-prop-dx}
```
  SReachTools/SimpleBox/dx  Box side half-lengths
  =================================================================
  
  Array of the half-lenghs of each side of the box
 
  ASCII Example of 2-d box, 'x' are vertices, 'c' is center
 
   x --------------------------------- x
   |                 |                 |
   |                                   |
   |                dx(2)              |
   |                                   |
   |                 |                 |
   |                                   |
   | ---- dx(1) ---- c ---- dx(1) ---- |
   |                                   |
   |                 |                 |
   |                                   |
   |                dx(2)              |
   |                                   |
   |                 |                 |
   |                                   |
   x --------------------------------- x
 
 
```

### Property: dim
{:#SimpleBox-prop-dim}
{% include important-note.html content="Property currently has no help documentation." %}

### Method: getPolyhedron
{:#SimpleBox-method-getPolyhedron}
```
  SReachTools/SimpleBox/getPolyhedron  Get Polyhedron form of box
  ================================================================
 
  Class method to get the MPT Polyhedron representation of the 
  SimpleBox object
 
  Usage
  -----
  simpbox = SimpleBox(0, [1, 1]);
  poly = simpbox.getPolyhedron();
 
  ================================================================
 
  poly = obj.getPolyhedron()
 
  Inputs: None
 
  Outputs:
  --------
  poly - Polyhedron object
 
  ================================================================
 
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
  
```

### Method: computeGaussianProbability
{:#SimpleBox-method-computeGaussianProbability}
```
  SReachTools/SimpleBox/computeGaussianProbability  Compute the likelihood
  for Gaussian to be in box
  =====================================================================
  
  Method to compute the likelihood for a Gaussian random variable to 
  lie in the given box (SimpleBox object)
 
  Usage:
  ------
  simpbox = SimpleBox(0, [1, 1]);
  p = simpbox.computeGaussianProbability(mvncdf(simpbox.vertices));
  
  =====================================================================
  
  p = obj.computeGaussianProbability(vertex_probabilities)
  
  Inputs:
  -------
    vertex_probabilities - Probabilities of Gaussian at each vertex
 
  Outputs:
  --------
    p - Probability
 
  =====================================================================
 
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
  
```

### Method: getBounds
{:#SimpleBox-method-getBounds}
```
  SimpleBox/getBounds  Get upper and lower bounds for simple box
  ====================================================================
 
  Get upper and lower bounds for SimpleBox object.
 
  Usage
  -----
  simpbox = SimpleBox(0, [1, 1]);
  [lb, ub] = simpbox.getBounds();
  
  ====================================================================
 
  [lb, ub] = obj.getBounds();
 
  Inputs: None
 
  Outputs:
  --------
  lb - Lower bounds
  ub - Upper bounds
 
  ====================================================================
 
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
  
```

