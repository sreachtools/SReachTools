---
layout: docs
title: Tube.m
---

<ul class="doc-list">
    <li class="doc-list"><a href="#Tube">Tube</a></li>
    <ul class="doc-list">
        <li><a href="#Tube-Tube">Constructor</a></li>
        <li>Properties</li>
        <ul class="doc-list">
            <li class="doc-list"><a href="#Tube-prop-dim">dim</a></li>
        </ul>
        <li>Methods</li>
        <ul class="doc-list">
            <li class="doc-list"><a href="#Tube-method-contains">contains</a></li>
            <li class="doc-list"><a href="#Tube-method-concat">concat</a></li>
            <li class="doc-list"><a href="#Tube-method-polyArray2Tube">polyArray2Tube</a></li>
        </ul>
    </ul>
</ul>

{:#Tube}
### Tube
```
  Create a tube object
  ==========================================================================
 
  Tube class
 
  Usage:
  ------
  % Three different calling mechanisms, reach-avoid problem
  tt = Tube('reach-avoid', safe_set, target_set, time_horizon);
  
  % Viability Problem
  tt = Tube('viability', safe_set, time_horizon);
 
  % Note that both of the above mechanisms will yield a tube of length
  % time_horizon+1 --- T_0, T_1, ..., T_{time_horizon}.
  
  % General tube
  % Can use general Polyhedron objects, for self-containment of the 
  % usage example will use empty Polyhedron objects
  tt = Tube(Polyhedron(), Polyhedron(), Polyhedron());
 
  Given a tube 'orig_tube' and a polytope 'int_polytope', you can also 
  construct a tube that is the result of the intersection of the polytope 
  and sets in the tube as
  tt = Tube('intersect', orig_tube, int_polytope);
    
  ==========================================================================
 
  Tube Properties:
  ------------------------
    dim - Dimension of the tube, i.e. dimension of the sets in the tube
 
  Tube Methods:
  ---------------------
    Tube/Tube       - Class constructor
    concat          - Get concatenated tube (Cartesian product of the polytopes)
    contains        - Check if a given concatenates state trajectory lies in the
                      tube
    polyArray2Tube  - Convert the array of polyhedra to a Tube object
  
  Apart from these, the following MATLAB functions have been overloaded
    disp            - Display the tube
    length          - Provide the length of the tube 
    end             - Indexing the end of the tube
    subsref         - Permits use of indexing of tube
 
  Notes:
  ------
  * MATLAB DEPENDENCY: MPT 3.0
  
  =========================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
  

    Reference page in Doc Center
       doc Tube

```

{:#Tube-Tube}
### Constructor
```
  Create a tube object
  ==========================================================================
 
  Tube class
 
  Usage:
  ------
  % Three different calling mechanisms, reach-avoid problem
  tt = Tube('reach-avoid', safe_set, target_set, time_horizon);
  
  % Viability Problem
  tt = Tube('viability', safe_set, time_horizon);
 
  % Note that both of the above mechanisms will yield a tube of length
  % time_horizon+1 --- T_0, T_1, ..., T_{time_horizon}.
  
  % General tube
  % Can use general Polyhedron objects, for self-containment of the 
  % usage example will use empty Polyhedron objects
  tt = Tube(Polyhedron(), Polyhedron(), Polyhedron());
 
  Given a tube 'orig_tube' and a polytope 'int_polytope', you can also 
  construct a tube that is the result of the intersection of the polytope 
  and sets in the tube as
  tt = Tube('intersect', orig_tube, int_polytope);
    
  ==========================================================================
 
  Tube Properties:
  ------------------------
    dim - Dimension of the tube, i.e. dimension of the sets in the tube
 
  Tube Methods:
  ---------------------
    Tube/Tube       - Class constructor
    concat          - Get concatenated tube (Cartesian product of the polytopes)
    contains        - Check if a given concatenates state trajectory lies in the
                      tube
    polyArray2Tube  - Convert the array of polyhedra to a Tube object
  
  Apart from these, the following MATLAB functions have been overloaded
    disp            - Display the tube
    length          - Provide the length of the tube 
    end             - Indexing the end of the tube
    subsref         - Permits use of indexing of tube
 
  Notes:
  ------
  * MATLAB DEPENDENCY: MPT 3.0
  
  =========================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
  

    Reference page in Doc Center
       doc Tube

```

### Property: dim
{:#Tube-prop-dim}
```
  Tube/dim
  ==================================================================
  
  Dimension of the tube / sets in the tube
  
```

### Method: contains
{:#Tube-method-contains}
```
  Check if a given concatenates state trajectory lies in the tube
  ======================================================================
 
  This method is a wrapper over MPT's Polyhedron/contains 
  
  Usage: See checkViaMonteCarloSims
 
  ======================================================================
 
  [contains_flag] = contains(obj,X);
  
  Inputs: 
  -------
 
    X             - Concatenated state vector (or collection of it
                    arranged columnwise)
     
  Outputs:
  --------
    contains_flag - Row vector of length equal to columns in X
 
  Notes:
  ------
  * This function is useful for Monte-Carlo simulation-based
    verification
  
  ======================================================================
  
  This function is part of the Stochastic Optimal Control Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
   
```

### Method: concat
{:#Tube-method-concat}
```
  Get concatenated tube (Cartesian product of the polytopes)
  ======================================================================
 
  This method computes the half-space representation of the concatenated 
  tube. When no arguments are provided, it returns 
  safe_set^{time_horizon} x target_set, a huge polyhedron in the
  (obj.dim x time_horizon)-dimensional Euclidean space.
 
  [concat_tube_A, concat_tube_b] = concat(obj);
 
  The output matrices satisfy the relation that the a concatenated 
  state vector X lies in the reach-avoid tube if and only if
  concat_tube_A * [initial_state;X] <= concat_tube_b 
 
  When arguments are specified, it provides the Cartesian of
  specific time splice mentioned. Specifically, if the half space
  representation of sets from t=3 to 5 is desired for a given
  tube of length 10 (sets are defined for t=0 to 9), we
  provide
 
  [concat_tube_A, concat_tube_b] = concat(obj, [4 6]);
  
  The +1 added is to account for MATLAB's indexing which begins
  from 1. In other words, provide the starting and ending index of
  interest with respect to the Tube array of polyhedrons.
 
  Usage: See getLowerBoundStochReachAvoid.
 
  ======================================================================
 
  [concat_tube_A, concat_tube_b] = concat(obj,varagin);
 
  Inputs:
  -------
    time_limits - A 1x2 vector [a b] with the 1 <= a,b <= length(obj).
                  If a>b, then empty matrices are returned.
     
  Outputs:
  --------
    concat_tube_A - State matrix concatenated for tube
    concat_tube_b - Input matrix concatenated for tube
 
  Notes:
  ------
  * This function also serves as a delegatee for input handling.
  
  ======================================================================
  
  This function is part of the Stochastic Optimal Control Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
 
```

### Method: polyArray2Tube
{:#Tube-method-polyArray2Tube}
```
  Convert the array of polyhedra to a Tube object
  ======================================================================
  Inputs:
  -------
    polyArray - Array of polyhedrons
 
  Outputs:
  --------
    newobj    - Tube object
 
  ======================================================================
  
  This function is part of the Stochastic Optimal Control Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
 
```

