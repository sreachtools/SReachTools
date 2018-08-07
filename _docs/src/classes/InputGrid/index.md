---
layout: docs
title: InputGrid.m
---

<ul class="doc-list">
    <li class="doc-list"><a href="#InputGrid">InputGrid</a></li>
    <ul class="doc-list">
        <li><a href="#InputGrid-InputGrid">Constructor</a></li>
        <li>Properties</li>
        <ul class="doc-list">
            <li class="doc-list"><a href="#InputGrid-prop-grid">grid</a></li>
            <li class="doc-list"><a href="#InputGrid-prop-lower_bounds">lower_bounds</a></li>
            <li class="doc-list"><a href="#InputGrid-prop-upper_bounds">upper_bounds</a></li>
            <li class="doc-list"><a href="#InputGrid-prop-n_points">n_points</a></li>
            <li class="doc-list"><a href="#InputGrid-prop-grid_delta">grid_delta</a></li>
            <li class="doc-list"><a href="#InputGrid-prop-dim">dim</a></li>
        </ul>
        <li>Methods</li>
        <ul class="doc-list">
            <li class="doc-list"><a href="#InputGrid-method-getIndicatorVectorForSet">getIndicatorVectorForSet</a></li>
        </ul>
    </ul>
</ul>

{:#InputGrid}
### InputGrid
```
  SReachTools/InputGrid: Create a input space grid object
  ============================================================================
 
  Defines an input space grid used for dynamic programming computations
 
  Usage:
  ------
 
  % Define a 2-dimensional state-space grid that extends from x = [-1, 1],
  % y = [-1, 1] with 100 points in each dimension
 
  grid = InputGrid([-1, -1], [1, 1], 100)
 
  % Can also define different dimensional spacings
 
  grid = InputGrid([-1, -1], [1, 1], [100, 50])
    
  ============================================================================
 
  InputGrid Properties:
  ---------------------
    grid         - Array of grid vectors, size prod(n_points) x dim
    lower_bounds - Lower bounds provided during construction
    upper_bounds - Upper bounds provided during construction
    n_points     - Number of points in grid in each dimension
    grid_delta   - Grid spacing, spacing between two grid points is 
                   2*grid_delta(i); 'i' being the dimension of interest
    dim          - Total number of dimensions in the grid
  
  InputGrid Methods:
  ------------------
    InputGrid/InputGrid      - Class constructor
    getIndicatorVectorForSet - Method to get an indicator vector of which grid 
                               points are in a Polyhedron set
 
  Notes:
  ------
  * After the 'external' option is removed from the SpaceGrid class there will
    be very little difference between InputGrid and SpaceGrid. The classes
    will be merged for convenience.
  
  ============================================================================
 
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  

    Reference page in Doc Center
       doc InputGrid

```

{:#InputGrid-InputGrid}
### Constructor
```
  SReachTools/InputGrid: Create a input space grid object
  ============================================================================
 
  Defines an input space grid used for dynamic programming computations
 
  Usage:
  ------
 
  % Define a 2-dimensional state-space grid that extends from x = [-1, 1],
  % y = [-1, 1] with 100 points in each dimension
 
  grid = InputGrid([-1, -1], [1, 1], 100)
 
  % Can also define different dimensional spacings
 
  grid = InputGrid([-1, -1], [1, 1], [100, 50])
    
  ============================================================================
 
  InputGrid Properties:
  ---------------------
    grid         - Array of grid vectors, size prod(n_points) x dim
    lower_bounds - Lower bounds provided during construction
    upper_bounds - Upper bounds provided during construction
    n_points     - Number of points in grid in each dimension
    grid_delta   - Grid spacing, spacing between two grid points is 
                   2*grid_delta(i); 'i' being the dimension of interest
    dim          - Total number of dimensions in the grid
  
  InputGrid Methods:
  ------------------
    InputGrid/InputGrid      - Class constructor
    getIndicatorVectorForSet - Method to get an indicator vector of which grid 
                               points are in a Polyhedron set
 
  Notes:
  ------
  * After the 'external' option is removed from the SpaceGrid class there will
    be very little difference between InputGrid and SpaceGrid. The classes
    will be merged for convenience.
  
  ============================================================================
 
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  

    Reference page in Doc Center
       doc InputGrid

```

### Property: grid
{:#InputGrid-prop-grid}
{% include important-note.html content="Property currently has no help documentation." %}

### Property: lower_bounds
{:#InputGrid-prop-lower_bounds}
{% include important-note.html content="Property currently has no help documentation." %}

### Property: upper_bounds
{:#InputGrid-prop-upper_bounds}
{% include important-note.html content="Property currently has no help documentation." %}

### Property: n_points
{:#InputGrid-prop-n_points}
{% include important-note.html content="Property currently has no help documentation." %}

### Property: grid_delta
{:#InputGrid-prop-grid_delta}
{% include important-note.html content="Property currently has no help documentation." %}

### Property: dim
{:#InputGrid-prop-dim}
{% include important-note.html content="Property currently has no help documentation." %}

### Method: getIndicatorVectorForSet
{:#InputGrid-method-getIndicatorVectorForSet}
```
  SReachTools/InputGrid/getIndicatorVectorForSet Get indicator vector
  for points in grid that lie in set 
  ====================================================================
 
  Class method returning an indicator vector (vector of zeros and ones)
  for all points in the grid that lie in a Polyhedral set
 
  Usage:
  ------
  grid = InputGrid([-1, -1], [1, 1], 100);
  s = Polyhedron('lb', [0, 0], 'ub', [1, 1]);
  vec = grid.getIndicatorVectorForSet(s);
  
  ====================================================================
 
  ind_vector = getIndicatorVectorForSet(obj, s)       
  
  Inputs:
  -------
    obj - InputGrid object
    s   - Polyhedron set
  
  Outputs:
  --------
    ind_vector - Indicator vector (n_points x 1)
  
  ====================================================================
 
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
```

