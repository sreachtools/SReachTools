---
layout: docs
title: SReachLagController.m
---

<ul class="doc-list">
    <li class="doc-list"><a href="#SReachLagController">SReachLagController</a></li>
    <ul class="doc-list">
        <li><a href="#SReachLagController-SReachLagController">Constructor</a></li>
        <li>Properties</li>
        <ul class="doc-list">
            <li class="doc-list"><a href="#SReachLagController-prop-tube">tube</a></li>
            <li class="doc-list"><a href="#SReachLagController-prop-system">system</a></li>
            <li class="doc-list"><a href="#SReachLagController-prop-dist_set">dist_set</a></li>
            <li class="doc-list"><a href="#SReachLagController-prop-time_horizon">time_horizon</a></li>
        </ul>
        <li>Methods</li>
        <ul class="doc-list">
            <li class="doc-list"><a href="#SReachLagController-method-getInput">getInput</a></li>
        </ul>
    </ul>
</ul>

{:#SReachLagController}
### SReachLagController
```
  Controller that maximizes stochastic reachability via Lagrangian-based
  computations and respects hard control bounds
  ==========================================================================
 
  SReachLagController class
 
  Usage:
  ------
  % Given a system sys, probability threshold prob_thresh, and a target_tube, we
  % can compute the Lagrangian-based underapproximation via these two commands.
  lagunder_options = SReachSetOptions('term', 'lag-under', ...
       'bound_set_method', 'ellipsoid', 'compute_style','vfmethod');
  [polytope_lagunder, extra_info_under] = SReachSet('term', 'lag-under', ...
          sys, prob_thresh, target_tube, lagunder_options);
  % We compute the associated controller using the following command
  srlcontrol = SReachLagController(sys, bounded_dist_set, robust_reach_tube);
  
  See also, example/cwhSReachSet.
 
  ==========================================================================
 
  SReachLagController Properties:
  -------------------------------
    system          - System under study [LtvSystem/LtiSystem object]
    dist_set        - Bounded disturbance set for which robustness computation
                      has been completed [Polyhedron/SReachEllipsoid object]
    tube            - Time-stamped robust sets in the state space from which
                      there is a control that works irrespective of the
                      disturbance in dist_set [Tube object]. See notes.
    time_horizon    - Computed from the obj.tube (Scalar)
 
  SReachLagController Methods:
  ----------------------------
    SReachLagController/SReachLagController 
                    - Constructor for SReachLagController
    SReachLagController/getInput 
                    - Get the input at time t given the state at time t, that
                      lies in the effective target tube
  
  Notes:
  ------
  * MATLAB DEPENDENCY: MPT 3.0
  * Robust target tube can also be provided as an array of polyhedron.
    The construct converts it into a tube object using
    Tube.polyArray2Tube().
  
  =========================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
  

    Reference page in Doc Center
       doc SReachLagController

```

{:#SReachLagController-SReachLagController}
### Constructor
```
  Controller that maximizes stochastic reachability via Lagrangian-based
  computations and respects hard control bounds
  ==========================================================================
 
  SReachLagController class
 
  Usage:
  ------
  % Given a system sys, probability threshold prob_thresh, and a target_tube, we
  % can compute the Lagrangian-based underapproximation via these two commands.
  lagunder_options = SReachSetOptions('term', 'lag-under', ...
       'bound_set_method', 'ellipsoid', 'compute_style','vfmethod');
  [polytope_lagunder, extra_info_under] = SReachSet('term', 'lag-under', ...
          sys, prob_thresh, target_tube, lagunder_options);
  % We compute the associated controller using the following command
  srlcontrol = SReachLagController(sys, bounded_dist_set, robust_reach_tube);
  
  See also, example/cwhSReachSet.
 
  ==========================================================================
 
  SReachLagController Properties:
  -------------------------------
    system          - System under study [LtvSystem/LtiSystem object]
    dist_set        - Bounded disturbance set for which robustness computation
                      has been completed [Polyhedron/SReachEllipsoid object]
    tube            - Time-stamped robust sets in the state space from which
                      there is a control that works irrespective of the
                      disturbance in dist_set [Tube object]. See notes.
    time_horizon    - Computed from the obj.tube (Scalar)
 
  SReachLagController Methods:
  ----------------------------
    SReachLagController/SReachLagController 
                    - Constructor for SReachLagController
    SReachLagController/getInput 
                    - Get the input at time t given the state at time t, that
                      lies in the effective target tube
  
  Notes:
  ------
  * MATLAB DEPENDENCY: MPT 3.0
  * Robust target tube can also be provided as an array of polyhedron.
    The construct converts it into a tube object using
    Tube.polyArray2Tube().
  
  =========================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
  

    Reference page in Doc Center
       doc SReachLagController

```

### Property: tube
{:#SReachLagController-prop-tube}
```
  SReachLagController/tube
  ==================================================================
  
  Time-stamped robust sets in the state space from which there is a
  control that works irrespective of the disturbance in dist_set [Tube
  object]. When Polyhedron array of sets for the tube is given, the
  constructor will use Tube.polyArray2Tube() to convert it into a tube.
  Half-space representation preferred.
 
```

### Property: system
{:#SReachLagController-prop-system}
```
  SReachLagController/system
  ==================================================================
  
  LtvSystem or LtiSystem object describing the dynamics
 
```

### Property: dist_set
{:#SReachLagController-prop-dist_set}
```
  SReachLagController/dist_set
  ==================================================================
  
  Polyhedron or SReachEllipsoid object describing the disturbance set
 
```

### Property: time_horizon
{:#SReachLagController-prop-time_horizon}
```
  SReachLagController/time_horizon
  ==================================================================
  
  Time horizon up to which the controller has been defined 
  (Defined using the length(obj.tube) and must be >0)
 
```

### Method: getInput
{:#SReachLagController-method-getInput}
```
  Get the input at time t given the state at time t, that lies in the
  effective target tube
  ======================================================================
  Inputs:
  -------
    obj             - SReachLagController object
    current_state   - Current state (a obj.system.state_dim x 1 vector)
    current_time    - Current time (a scalar integer) 
 
  Outputs:
  --------
    action          - Action to apply at time t
 
  Notes:
  ------
  * If infeasible initial state is provided or the feasible input space
    turns out to be empty, an Invalid Arguments error is thrown.
    Therefore, while using this function to generate trajectories, make
    sure that the disturbances lie in obj.dist_set.
 
  ======================================================================
  
  This function is part of the Stochastic Optimal Control Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
 
```

