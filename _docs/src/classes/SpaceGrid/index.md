---
layout: docs
title: SpaceGrid.m
---

<ul class="doc-list">
    <li class="doc-list"><a href="#SpaceGrid">SpaceGrid</a></li>
    <ul class="doc-list">
        <li><a href="#SpaceGrid-SpaceGrid">Constructor</a></li>
        <li>Properties</li>
        <ul class="doc-list">
            <li class="doc-list"><a href="#SpaceGrid-prop-grid">grid</a></li>
            <li class="doc-list"><a href="#SpaceGrid-prop-lower_bounds">lower_bounds</a></li>
            <li class="doc-list"><a href="#SpaceGrid-prop-upper_bounds">upper_bounds</a></li>
            <li class="doc-list"><a href="#SpaceGrid-prop-n_points">n_points</a></li>
            <li class="doc-list"><a href="#SpaceGrid-prop-grid_delta">grid_delta</a></li>
            <li class="doc-list"><a href="#SpaceGrid-prop-dim">dim</a></li>
        </ul>
        <li>Methods</li>
        <ul class="doc-list">
            <li class="doc-list"><a href="#SpaceGrid-method-plotGridProbability">plotGridProbability</a></li>
            <li class="doc-list"><a href="#SpaceGrid-method-getExternalGrid">getExternalGrid</a></li>
            <li class="doc-list"><a href="#SpaceGrid-method-getMeshGrids">getMeshGrids</a></li>
            <li class="doc-list"><a href="#SpaceGrid-method-sortGrid">sortGrid</a></li>
            <li class="doc-list"><a href="#SpaceGrid-method-getIndicatorVectorForSet">getIndicatorVectorForSet</a></li>
        </ul>
    </ul>
</ul>

{:#SpaceGrid}
### SpaceGrid
```
  SReachTools/SpaceGrid  Create a state space grid object
  =============================================================================
 
  Class to hold the gridding of a particular space, e.g. state or input.
 
  Usage:
  ------
  % Define a 2-dimensional state-space grid that extends from x = [-1, 1],
  % y = [-1, 1] with 100 points in each dimension
 
  grid = SpaceGrid([-1, -1], [1, 1], 100)
 
  % Can also define different dimensional spacings
 
  grid = SpaceGrid([-1, -1], [1, 1], [100, 50])
  
  =============================================================================
 
  SpaceGrid Properties:
  ---------------------
    grid         - Array of grid vectors, size prod(n_points) x dim
    lower_bounds - Lower bounds provided during construction
    upper_bounds - Upper bounds provided during construction
    n_points     - Number of points in grid in each dimension
    grid_delta   - Grid spacing, spacing between two grid points is 
                   2*grid_delta(i); 'i' being the dimension of interest
    dim          - Total number of dimensions in the grid
  
  SpaceGrid Methods:
  ------------------
    SpaceGrid/SpaceGrid      - Constructor
    getIndicatorVectorForSet - Method to get an indicator vector of which grid 
                               points are in a Polyhedron set
    getMeshGrids             - Get the associated MATLAB mesh grids for the 
                               space; only works for 2 or 3-dimensional systems
    getExternalGrid          - Method to get an external grid for the current
                               grid
    plotGridProbability      - Helper method for plotting dynamic programming
                               probabilities on the grid
  
  =============================================================================
 
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
  

    Reference page in Doc Center
       doc SpaceGrid

```

{:#SpaceGrid-SpaceGrid}
### Constructor
```
  SReachTools/SpaceGrid  Create a state space grid object
  =============================================================================
 
  Class to hold the gridding of a particular space, e.g. state or input.
 
  Usage:
  ------
  % Define a 2-dimensional state-space grid that extends from x = [-1, 1],
  % y = [-1, 1] with 100 points in each dimension
 
  grid = SpaceGrid([-1, -1], [1, 1], 100)
 
  % Can also define different dimensional spacings
 
  grid = SpaceGrid([-1, -1], [1, 1], [100, 50])
  
  =============================================================================
 
  SpaceGrid Properties:
  ---------------------
    grid         - Array of grid vectors, size prod(n_points) x dim
    lower_bounds - Lower bounds provided during construction
    upper_bounds - Upper bounds provided during construction
    n_points     - Number of points in grid in each dimension
    grid_delta   - Grid spacing, spacing between two grid points is 
                   2*grid_delta(i); 'i' being the dimension of interest
    dim          - Total number of dimensions in the grid
  
  SpaceGrid Methods:
  ------------------
    SpaceGrid/SpaceGrid      - Constructor
    getIndicatorVectorForSet - Method to get an indicator vector of which grid 
                               points are in a Polyhedron set
    getMeshGrids             - Get the associated MATLAB mesh grids for the 
                               space; only works for 2 or 3-dimensional systems
    getExternalGrid          - Method to get an external grid for the current
                               grid
    plotGridProbability      - Helper method for plotting dynamic programming
                               probabilities on the grid
  
  =============================================================================
 
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
  

    Reference page in Doc Center
       doc SpaceGrid

```

### Property: grid
{:#SpaceGrid-prop-grid}
{% include important-note.html content="Property currently has no help documentation." %}

### Property: lower_bounds
{:#SpaceGrid-prop-lower_bounds}
{% include important-note.html content="Property currently has no help documentation." %}

### Property: upper_bounds
{:#SpaceGrid-prop-upper_bounds}
{% include important-note.html content="Property currently has no help documentation." %}

### Property: n_points
{:#SpaceGrid-prop-n_points}
{% include important-note.html content="Property currently has no help documentation." %}

### Property: grid_delta
{:#SpaceGrid-prop-grid_delta}
{% include important-note.html content="Property currently has no help documentation." %}

### Property: dim
{:#SpaceGrid-prop-dim}
{% include important-note.html content="Property currently has no help documentation." %}

### Method: plotGridProbability
{:#SpaceGrid-method-plotGridProbability}
```
  SReachTools/SpaceGrid/plotGridProbability  Plot grid probability
  ====================================================================
  
  Perform surface plot of 2-dimensional grid probability
  
  Usage:
    grid = SpaceGrid([-1, -1], [1, 1], 100);
    grid.plotGridProbability(mvncdf(grid.grid));
  
  ====================================================================
  
  obj.plotGridProbability
  
  Inputs:  None
  Outputs: None
  
  ====================================================================
 
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```

### Method: getExternalGrid
{:#SpaceGrid-method-getExternalGrid}
{% include important-note.html content="Method currently has no help documentation." %}

### Method: getMeshGrids
{:#SpaceGrid-method-getMeshGrids}
```
  SReachTools/SpaceGrid/getMeshGrids  Get MATLAB meshgrids
  ====================================================================
  
  Get MATLAB meshgrids for the SpaceGrid object; only works for grids 
  that are 2 or 3-dimensional
  
  Usage:
  ------
    grid = SpaceGrid([-1, -1], [1, 1], 100);
    [X,Y] = grid.getMeshGrids();
  
  ====================================================================
  
  [X,Y] = obj.getMeshGrids();
  [X,Y,Z] = obj.getMeshGrids();
 
  Inputs: None
  
  Outputs:
  --------
    X - x-meshgrid
    Y - y-meshgrid
    Z - z-meshgrid
 
  ====================================================================
 
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```

### Method: sortGrid
{:#SpaceGrid-method-sortGrid}
```
  SReachTools/SpaceGrid/sortGrid  Sort space grid vectors
  ====================================================================
  
  Sort space grid vectors, ascending
 
  Usage:
  ------
    grid = SpaceGrid([-1, -1], [1, 1], 100);
    grid.sortGrid();
 
  ====================================================================
 
  sortGrid(obj)
 
  Inputs:  None
  Outputs: None
  
  ====================================================================
 
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```

### Method: getIndicatorVectorForSet
{:#SpaceGrid-method-getIndicatorVectorForSet}
```
  SReachTools/SpaceGrid/getIndicatorVectorForSet  Get indicator vector
  ====================================================================
  
  Get indicator vector for the grid points which lie in a Polyhedron
  set, s.
 
  Usage:
  ------
    grid = SpaceGrid([-1, -1], [1, 1], 100);
    ind_vector = grid.getIndicatorVectorForSet(Polyhedron(...
        'lb', [0,0], 'ub', [1, 1]));
 
  ====================================================================
  
  ind_vector = getIndicatorVectorForSet(obj, s)
 
  Inputs:
  -------
    s - Polyhedron object set
  
  Outputs:
  --------
    ind_vector - Indicator vector for grid points in set s
 
  ====================================================================
 
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
  
```

