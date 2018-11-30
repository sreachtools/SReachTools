---
layout: docs
title: getDubinsCarLtv.m
---

```
  Get a LtvSystem object for the Dubins car model with known turning rate
  sequence
  ============================================================================
  
  Usage:
  ------
  % Known turning rate sequence
  time_horizon = 50;
  omega = pi/time_horizon/sampling_time;
  turning_rate = omega*ones(time_horizon,1);
  init_heading = pi/10;                       % Initial heading 
  sampling_time = 0.1;                        % Sampling time
  % Input space definition
  umax = 6;
  input_space = Polyhedron('lb',0,'ub',umax);
  % Disturbance matrix and random vector definition
  dist_matrix = eye(2);
  eta_dist = RandomVector('Gaussian',zeros(2,1), 0.001 * eye(2));
  
  [sys, heading_vec] = getDubinsCarLtv('add-dist', turning_rate, ...
    init_heading, sampling_time, input_space, dist_matrix, eta_dist);
 
  ============================================================================
  
  sys = getDubinsCarLtv(type, turning_rate_seq, initial_heading, sampling_time,
    varargin)
 
  Inputs:
  -------
    type        - 
    turning_rate_seq    - Known turning rate sequence (column vector of length
                          time_horizon) 
    initial_heading     - Initial heading angle
    sampling_time       - Sampling time for the system
    Required additional arguments for different types:
    type: 'add-dist' (Additive disturbance with velocity as input)
        velocity_input  - Bounds on the velocity (1-dimensional Polyhedron
                          object)
        dist_matrix     - Disturbance matrix for the additive disturbance 
        dist            - Disturbance
    type = 'vel-dist' (Velocity is the disturbance) 
        velocity_dist   - Velocity disturbance
  
  Outputs:
  --------
    sys                 - LtvSystem object describing the Dubin's car
  
  =============================================================================
 
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/abyvinod/SReachTools/blob/master/LICENSE
  
 
```
