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
            <li class="doc-list"><a href="#LtvSystem-method-getHmatMeanCovForXSansInput">getHmatMeanCovForXSansInput</a></li>
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
  SReachTools/LtvSystem: Create a discrete-time LTV system object
  ============================================================================
 
  Defines a discrete-time LTV system that is:
      - control-free and disturbance-free, or
      - controlled but disturbance-free, or
      - perturbed (stochastic/uncertain) but control-free, or
      - controlled and perturbed (stochastic/uncertain).
 
  Perturbation can be either:
      - a bounded uncertainity with no stochastic information
      - a StochasticDisturbance object
 
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
    getHmatMeanCovForXSansInput
                          - Get input policy-free mean and covariance of the
                            trajectory from a given initial state for a known
                            time horizon and the concatenated input matrix
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
  SReachTools/LtvSystem: Create a discrete-time LTV system object
  ============================================================================
 
  Defines a discrete-time LTV system that is:
      - control-free and disturbance-free, or
      - controlled but disturbance-free, or
      - perturbed (stochastic/uncertain) but control-free, or
      - controlled and perturbed (stochastic/uncertain).
 
  Perturbation can be either:
      - a bounded uncertainity with no stochastic information
      - a StochasticDisturbance object
 
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
    getHmatMeanCovForXSansInput
                          - Get input policy-free mean and covariance of the
                            trajectory from a given initial state for a known
                            time horizon and the concatenated input matrix
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
{% include important-note.html content="Property currently has no help documentation." %}

### Property: state_dim
{:#LtvSystem-prop-state_dim}
{% include important-note.html content="Property currently has no help documentation." %}

### Property: input_mat
{:#LtvSystem-prop-input_mat}
{% include important-note.html content="Property currently has no help documentation." %}

### Property: input_space
{:#LtvSystem-prop-input_space}
{% include important-note.html content="Property currently has no help documentation." %}

### Property: input_dim
{:#LtvSystem-prop-input_dim}
{% include important-note.html content="Property currently has no help documentation." %}

### Property: dist
{:#LtvSystem-prop-dist}
{% include important-note.html content="Property currently has no help documentation." %}

### Property: dist_mat
{:#LtvSystem-prop-dist_mat}
{% include important-note.html content="Property currently has no help documentation." %}

### Property: dist_dim
{:#LtvSystem-prop-dist_dim}
{% include important-note.html content="Property currently has no help documentation." %}

### Method: getHmatMeanCovForXSansInput
{:#LtvSystem-method-getHmatMeanCovForXSansInput}
```
  SReachTools/LtvSystem/getHmatMeanCovForXSansInput: Get input policy-free mean
  and covariance of the trajectory from a given initial state for a known time
  horizon and the concatenated input matrix
  ============================================================================
 
  Helps in the computation of the mean and covariance of the concatenated
  state vector X for a given stochastic LTI system as given in (17) of
 
  A. Vinod and M. Oishi, "Scalable Underapproximation for Stochastic Reach-Avoid
  Problem for High-Dimensional LTI Systems using Fourier Transforms," in IEEE
  Control Systems Letters (L-CSS), 2017.
 
  Also, returns H, and Z and G if needed
 
  For more details on the matrix notation, please see the documentation of
  LtvSystem/getConcatMats(). 
 
  Usage: See getLowerBoundStochReachAvoid
 
  ============================================================================
  
  [H, mean_X_sans_input, cov_X_sans_input, varargout] = ...
                 getHmatMeanCovForXSansInput(sys, ...
                                             initial_state, ...
                                             time_horizon)
  Inputs:
  -------
    sys           - An object of LtvSystem class 
    initial_state - Initial state can be a deterministic n-dimensional vector
                    x_0 or a RandomVector object
    time_horizon  - Time of interest (N)
 
  Outputs:
  --------
    H                - Concatenated input matrix
    mean_X_sans_input- Mean of X with zero input under the disturbance from the
                       provided initial state
    cov_X_sans_input - Covariance of X with zero input under the disturbance
                       from the provided initial state
    Z                - (optional) Concatenated state matrix
    G                - (optional) Concatenated disturbance matrix
 
  Notes:
  ------
  * X refers to the concatenated state vector X=[x_1^\top x_2^\top ...
    x_N^\top]^\top. See @LtvSystem/getConcatMats for more
    information about the notation used.
  * Assumes the disturbance is independent and identically distributed
  * This function also serves as a delegatee for input handling.
  
  ============================================================================
 
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
 
 
```

### Method: getConcatInputSpace
{:#LtvSystem-method-getConcatInputSpace}
```
  SReachTools/LtvSystem/getConcatInputSpace: Get half space representation of
  the concatenated (polytopic) input space for the given time horizon
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
  SReachTools/LtvSystem/getConcatMats: Get concatenated matrices
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
  % Compute the robust reach-avoid set
  [Z,H,G] = getConcatMats(sys, time_horizon);
 
  =============================================================================
 
  [Z,H,G] = getConcatMats(sys, time_horizon)
  Inputs:
  -------
    sys          - An object of LtvSystem class 
    time_horizon - Time of interest (N)
 
  Outputs:
  --------
    Z - Concatenated state matrix
    H - Concatenated input matrix
    G - Concatenated disturbance matrix
 
  Notes:
  ------
  * For control-free and/or disturbance-free LTI systems, H and G are set to
    zeros( sys.state_dim * time_horizon, 1) as appropriate.
  * Deviation from Skaf and Boyd's definition,
      * Concatenated state is X=[x_1 x_2 ... x_{N}].
      * Z definition excludes the initial state x_0 in contrast to Skaf and
        Boyd's definition of x_0.
      * H, G does include the first row since initial state is not there.
      * This function computes for a LTI system instead of the original LTV
        formulation.
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
  SReachTools/LtvSystem/islti: Get boolean result if system is LTI
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
  SReachTools/LtvSystem/islti: Get boolean result if system is LTI
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

