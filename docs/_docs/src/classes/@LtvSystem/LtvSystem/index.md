---
layout: docs
title: LtvSystem.m
---

<ul class="doc-list">
    <li class="doc-list"><a href="#LtvSystem">LtvSystem</a></li>
    <ul class="doc-list">
        <li><a href="#LtvSystem-LtvSystem">Constructor</a></li>
        <li>Properties</li>
        <ul class="doc-list">
            <li class="doc-list"><a href="#LtvSystem-prop-state_mat">state_mat</a></li>
            <li class="doc-list"><a href="#LtvSystem-prop-state_dim">state_dim</a></li>
            <li class="doc-list"><a href="#LtvSystem-prop-input_mat">input_mat</a></li>
            <li class="doc-list"><a href="#LtvSystem-prop-input_space">input_space</a></li>
            <li class="doc-list"><a href="#LtvSystem-prop-input_dim">input_dim</a></li>
            <li class="doc-list"><a href="#LtvSystem-prop-dist">dist</a></li>
            <li class="doc-list"><a href="#LtvSystem-prop-dist_mat">dist_mat</a></li>
            <li class="doc-list"><a href="#LtvSystem-prop-dist_dim">dist_dim</a></li>
        </ul>
        <li>Methods</li>
        <ul class="doc-list">
            <li class="doc-list"><a href="#LtvSystem-method-getConcatInputSpace">getConcatInputSpace</a></li>
            <li class="doc-list"><a href="#LtvSystem-method-getConcatMats">getConcatMats</a></li>
            <li class="doc-list"><a href="#LtvSystem-method-isltv">isltv</a></li>
            <li class="doc-list"><a href="#LtvSystem-method-islti">islti</a></li>
        </ul>
    </ul>
</ul>

{:#LtvSystem}
### LtvSystem
```
  Create a discrete-time LTV system object
  ============================================================================
 
  Defines a discrete-time LTV system that is:
      - control-free and disturbance-free, or
      - controlled but disturbance-free, or
      - perturbed (stochastic/uncertain) but control-free, or
      - controlled and perturbed (stochastic/uncertain).
 
  Perturbation can be either:
      - a bounded uncertainity with no stochastic information
      - a RandomVector object
 
   Usage:
   ------
   % Define a double integrator system:
 
   T = 0.5;
   sys = LtvSystem('StateMatrix', [1, T; 0, 1], ...
                   'InputMatrix', [T^2/2;T], ...
                   'InputSpace', Polyhedron('lb', -1, 'ub', 1), ...
                   'DisturbanceMatrix', [T^2/2;T], ...
                   'Disturbance', Polyhedron('lb', -1, 'ub', 1));
    
  ============================================================================
 
  LtvSystem Properties:
  ---------------------
    state_mat       - State matrix (Square matrix, state_dim x state_dim)
    input_mat       - Input matrix (Matrix, state_dim x input_dim)
    input_space     - Input space (empty / Polyhedron)
    dist_mat        - Disturbance matrix (Matrix, state_dim x dist_dim)
    dist            - Disturbance object (empty/Polyhedron/RandomVector)     
    state_dim       - State dimension (scalar)   
    input_dim       - Input dimension (scalar)  
    dist_dim        - Disturbance dimension (scalar)
  
  LtvSystem Methods:
  ------------------
    LtvSystem/LtvSystem   - Constructor
    getConcatInputSpace   - Get concatenated input space
    getConcatMats         - Get concatenated state, input, and disturbance
                            matrices
    islti                 - Get logical value 1 if system is LTI
    isltv                 - Get logical value 1 if system is LTV (strictly)
  
  Notes:
  ------
  * EXTERNAL DEPENDENCY: Uses MPT3 to define input,robust disturbance space
 
  =============================================================================
 
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
 

    Reference page in Doc Center
       doc LtvSystem

```

{:#LtvSystem-LtvSystem}
### Constructor
```
  Create a discrete-time LTV system object
  ============================================================================
 
  Defines a discrete-time LTV system that is:
      - control-free and disturbance-free, or
      - controlled but disturbance-free, or
      - perturbed (stochastic/uncertain) but control-free, or
      - controlled and perturbed (stochastic/uncertain).
 
  Perturbation can be either:
      - a bounded uncertainity with no stochastic information
      - a RandomVector object
 
   Usage:
   ------
   % Define a double integrator system:
 
   T = 0.5;
   sys = LtvSystem('StateMatrix', [1, T; 0, 1], ...
                   'InputMatrix', [T^2/2;T], ...
                   'InputSpace', Polyhedron('lb', -1, 'ub', 1), ...
                   'DisturbanceMatrix', [T^2/2;T], ...
                   'Disturbance', Polyhedron('lb', -1, 'ub', 1));
    
  ============================================================================
 
  LtvSystem Properties:
  ---------------------
    state_mat       - State matrix (Square matrix, state_dim x state_dim)
    input_mat       - Input matrix (Matrix, state_dim x input_dim)
    input_space     - Input space (empty / Polyhedron)
    dist_mat        - Disturbance matrix (Matrix, state_dim x dist_dim)
    dist            - Disturbance object (empty/Polyhedron/RandomVector)     
    state_dim       - State dimension (scalar)   
    input_dim       - Input dimension (scalar)  
    dist_dim        - Disturbance dimension (scalar)
  
  LtvSystem Methods:
  ------------------
    LtvSystem/LtvSystem   - Constructor
    getConcatInputSpace   - Get concatenated input space
    getConcatMats         - Get concatenated state, input, and disturbance
                            matrices
    islti                 - Get logical value 1 if system is LTI
    isltv                 - Get logical value 1 if system is LTV (strictly)
  
  Notes:
  ------
  * EXTERNAL DEPENDENCY: Uses MPT3 to define input,robust disturbance space
 
  =============================================================================
 
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
 

    Reference page in Doc Center
       doc LtvSystem

```

### Property: state_mat
{:#LtvSystem-prop-state_mat}
```
  LTVSystem/state_mat
  ====================================================================
  
  System state matrix
    - MATLAB matrix for LTI system
    - Anonymous functions for LTV system, maps time index to matrix
  
```

### Property: state_dim
{:#LtvSystem-prop-state_dim}
```
  LTVSystem/state_dim
  ====================================================================
  
  Dimension of the state vector
  
```

### Property: input_mat
{:#LtvSystem-prop-input_mat}
```
  LTVSystem/input_mat
  ====================================================================
  
  System input matrix
    - MATLAB matrix for LTI system
    - Anonymous functions for LTV system, maps time index to matrix
  
```

### Property: input_space
{:#LtvSystem-prop-input_space}
```
  LTVSystem/input_space
  ====================================================================
  
  System input space, polyhedron object
  
```

### Property: input_dim
{:#LtvSystem-prop-input_dim}
```
  LTVSystem/input_dim
  ====================================================================
  
  Dimension of the input vector
  
```

### Property: dist
{:#LtvSystem-prop-dist}
```
  LTVSystem/dist
  ====================================================================
  
  LTV system disturbance, RandomVector object
  
```

### Property: dist_mat
{:#LtvSystem-prop-dist_mat}
```
  LTVSystem/dist_mat
  ====================================================================
  
  System disturbance matrix
    - MATLAB matrix for LTI system
    - Anonymous functions for LTV system, maps time index to matrix
  
```

### Property: dist_dim
{:#LtvSystem-prop-dist_dim}
```
  LTVSystem/dist_dim
  ====================================================================
  
  Dimension of disturbance vector
  
```

### Method: getConcatInputSpace
{:#LtvSystem-method-getConcatInputSpace}
```
  Get half space representation of the concatenated (polytopic) input space 
  for the given time horizon
  ============================================================================
  
  Computes the input_space^{time_horizon} corresponding to a given set, which
  is the set of admissible open-loop control polices. This function computes the
  half-space representation of the cartesian products of polytopic input spaces.
 
  Usage:
  ------
 
  % Compute the (matrix form) set of admissible open-loop control policies given
  % a LtvSystem and a time horizon
 
  sys = LtvSystem(...
      'StateMatrix', eye(2), ...
      'InputMatrix', ones(2,1), ...
      'InputSpace', Polyhedron('lb', -umax, 'ub', umax));
  time_horizon = 10;
  [concat_input_space_A, concat_input_space_b] = ...
                                              getConcatInputSpace(sys, ...
                                                                  time_horizon);
  
  ============================================================================
 
  [concat_input_space_A, concat_input_space_b] =...
                                               getConcatInputSpace(sys, ...
                                                                   time_horizon)
  
  Inputs:
  -------
    sys                  - An object of LtvSystem class 
    time_horizon         - Time horizon
 
  Outputs:
  --------
    concat_input_space_A, concat_input_space_b 
                         - Concatenated input space (Halfspace representation)
 
  =============================================================================
 
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```

### Method: getConcatMats
{:#LtvSystem-method-getConcatMats}
```
  Get concatenated matrices
  ============================================================================
  
  Computes the matrices corresponding to the concatentated state vector X.
 
  Consider a LtvSystem object with n as the state_dim, m as the
  input_dim, and p as the disturbance_dim. Given a time of
  interest N, we define a concatenated state vector (a nN-dimensional vector)
            __       __
            |   x_1   |
            |   x_2   |
        X = |   ...   |
            | x_{N-1} |          
            |  x_{N}  |          
            ---     ---
  where x_t is the state of the system with 1 <= t <= N.  Similarly, one can
  define concated input and noise vectors U and W (mN-dimensional and
  pN-dimensional vectors),
            __       __         __       __
            |   u_0   |         |   w_0   |
            |   u_1   |         |   w_1   |
        U = |   ...   |,   W  = |   ...   |
            | u_{N-2} |         | w_{N-2} |      
            | u_{N-1} |         | w_{N-1} |      
            ---     ---         ---     ---
 
  Given the initial state x_0, we have
 
        X = Z * x_0 + H * U + G * W
 
  where Z (nN x n matrix), H (nN x mN matrix), and G  (nN x
  pN matrix) are appropriate matrices. These matrices (with minor
  modifications noted below) are given in (3) in 
     J. Skaf and S. Boyd, "Design of Affine Controllers via Convex
     Optimization", in IEEE Trans. Automatic Control, 2010. 
 
  This function computes Z, H, and G.
 
  Usage:
  ------
 
  % Compute the concatenated matrices for a double integrator with a time of
  % interest, 10
 
  % Problem parameters
  time_horizon = 10;
  T = 0.25;
  umax = 0.75;
  dmax = 0.1;
  % Double integrator system
  sys = LtvSystem(...
      'StateMatrix', [1, T; 0, 1], ...
      'InputMatrix', [T^2; T], ...
      'InputSpace', Polyhedron('lb', -umax, 'ub', umax), ...
      'DisturbanceMatrix', eye(2), ...
      'Disturbance', Polyhedron('lb', -dmax *ones(2,1), 'ub', dmax *ones(2,1)));
  
  % Get the concatenated matrices
  [Z,H,G] = getConcatMats(sys, time_horizon);
 
  =============================================================================
 
  [Z,H,G] = getConcatMats(sys, time_horizon)
 
  Inputs:
  -------
    sys          - An object of LtvSystem class 
    time_horizon - Time horizon (N) with the control provided from 0 to N-1
 
  Outputs:
  --------
    Z - Concatenated state matrix
    H - Concatenated input matrix
    G - Concatenated disturbance matrix
 
  Notes:
  ------
  * For control-free and/or disturbance-free LTI/LTV systems, H and G are set to
    zeros( sys.state_dim * time_horizon, 1) as appropriate.
  * Deviation from Skaf and Boyd's definition,
      * Concatenated state is X=[x_1 x_2 ... x_{N}], with the initial state
      x_0 EXCLUDED in contrast to Skaf and Boyd's TAC 2010 definitions.
  * Computes the extended controllability matrix via for loops. (suboptimal way)
  * This function also serves as a delegatee for input handling
  
  ============================================================================
 
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```

### Method: isltv
{:#LtvSystem-method-isltv}
```
   Get boolean result if system is LTI
  ====================================================================
 
  Get boolean result if system is LTI. Considered LTI if state_mat, 
  input_mat, and dist_mat are all matrices
 
  Usage:
  ------
  T = 0.5;
  sys = LtvSystem('StateMatrix', [1, T; 0, 1], ...
                  'InputMatrix', [T^2/2;T], ...
                  'InputSpace', Polyhedron('lb', -1, 'ub', 1), ...
                  'DisturbanceMatrix', [T^2/2;T], ...
                  'Disturbance', Polyhedron('lb', -1, 'ub', 1));
  if sys.islti()
    disp('System is LTI')
  else
    disp('System is not LTI')
  end
  
  =====================================================================
 
  yn = obj.islti()
  
  Inputs: None
  Outputs:
  --------
    yn - Logical value of 1 if system is LTI
  
  =====================================================================
  
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
```

### Method: islti
{:#LtvSystem-method-islti}
```
  Get boolean result if system is LTI
  ====================================================================
 
  Get boolean result if system is LTI. Considered LTI if state_mat, 
  input_mat, and dist_mat are all matrices
 
  Usage:
  ------
  T = 0.5;
  sys = LtvSystem('StateMatrix', [1, T; 0, 1], ...
                  'InputMatrix', [T^2/2;T], ...
                  'InputSpace', Polyhedron('lb', -1, 'ub', 1), ...
                  'DisturbanceMatrix', [T^2/2;T], ...
                  'Disturbance', Polyhedron('lb', -1, 'ub', 1));
  if sys.islti()
    disp('System is LTI')
  else
    disp('System is not LTI')
  end
  
  =====================================================================
 
  yn = obj.islti()
  
  Inputs: None
  Outputs:
  --------
    yn - Logical value of 1 if system is LTI
  
  =====================================================================
  
    This function is part of the Stochastic Reachability Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
```

