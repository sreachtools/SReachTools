---
layout: docs
title: TargetTube.m
---

<ul class="doc-list">
    <li class="doc-list"><a href="#TargetTube">TargetTube</a></li>
    <ul class="doc-list">
        <li><a href="#TargetTube-TargetTube">Constructor</a></li>
        <li>Properties</li>
        <ul class="doc-list">
            <li class="doc-list"><a href="#TargetTube-prop-dim">dim</a></li>
        </ul>
        <li>Methods</li>
        <ul class="doc-list">
            <li class="doc-list"><a href="#TargetTube-method-contains">contains</a></li>
            <li class="doc-list"><a href="#TargetTube-method-concat">concat</a></li>
        </ul>
    </ul>
</ul>

{:#TargetTube}
### TargetTube
```
  SReachTools/TargetTube: Create a target tube object
  ==========================================================================
 
  Target tube class
 
  Usage:
  ------
  % Three different calling mechanisms, reach-avoid problem
  tt = TargetTube('reach-avoid', safe_set, target_set, time_horizon);
  
  % Viability Problem
  tt = TargetTube('viability', safe_set, time_horizon);
 
  Note that both of the above mechanisms will yield a target tube of length
  time_horizon+1 --- T_0, T_1, ..., T_{time_horizon}.
  
  % General tube
  % Can use general Polyhedron objects, for self-containment of the 
  % usage example will use empty Polyhedron objects
  tt = TargetTube(Polyhedron(), Polyhedron(), Polyhedron());
    
  ==========================================================================
 
  TargetTube Properties:
  ------------------------
 
  TargetTube Methods:
  ---------------------
    TargetTube/TargetTube - Class constructor
  
  Notes:
  ------
  * MATLAB DEPENDENCY: MPT 3.0
  
  =========================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
  

    Reference page in Doc Center
       doc TargetTube

```

{:#TargetTube-TargetTube}
### Constructor
```
  SReachTools/TargetTube: Create a target tube object
  ==========================================================================
 
  Target tube class
 
  Usage:
  ------
  % Three different calling mechanisms, reach-avoid problem
  tt = TargetTube('reach-avoid', safe_set, target_set, time_horizon);
  
  % Viability Problem
  tt = TargetTube('viability', safe_set, time_horizon);
 
  Note that both of the above mechanisms will yield a target tube of length
  time_horizon+1 --- T_0, T_1, ..., T_{time_horizon}.
  
  % General tube
  % Can use general Polyhedron objects, for self-containment of the 
  % usage example will use empty Polyhedron objects
  tt = TargetTube(Polyhedron(), Polyhedron(), Polyhedron());
    
  ==========================================================================
 
  TargetTube Properties:
  ------------------------
 
  TargetTube Methods:
  ---------------------
    TargetTube/TargetTube - Class constructor
  
  Notes:
  ------
  * MATLAB DEPENDENCY: MPT 3.0
  
  =========================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
  

    Reference page in Doc Center
       doc TargetTube

```

### Property: dim
{:#TargetTube-prop-dim}
{% include important-note.html content="Property currently has no help documentation." %}

### Method: contains
{:#TargetTube-method-contains}
```
  SReachTools/TargetTube/contains: Check if a given concatenates state
  trajectory lies in the target tube
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
{:#TargetTube-method-concat}
```
  SReachTools/TargetTube/concat: Get concatenated target tube
  ======================================================================
 
  This method computes the concatenated target tube, 
  safe_set^{time_horizon -1 } x target_set, a huge polyhedron in the
  (obj.dim x time_horizon)-dimensional Euclidean space.
  The output matrices satisfy the relation that the a concatenated 
  state vector X lies in the reach-avoid tube if and only if
  
  concat_target_tube_A * X <= concat_target_tube_b 
 
  Usage: See getFtLowerBoundTargetTube.
 
  ======================================================================
 
  [concat_target_tube_A, concat_target_tube_b] = concat(obj);
  
  Inputs: None
     
  Outputs:
  --------
    concat_target_tube_A - State matrix concatenated for target tube
    concat_target_tube_b - Input matrix concatenated for target tube
 
  Notes:
  ------
  * This function also serves as a delegatee for input handling.
  
  ======================================================================
  
  This function is part of the Stochastic Optimal Control Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
 
```

