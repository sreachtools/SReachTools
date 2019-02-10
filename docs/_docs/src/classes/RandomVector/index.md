---
layout: docs
title: RandomVector.m
---

<ul class="doc-list">
    <li class="doc-list"><a href="#RandomVector">RandomVector</a></li>
    <ul class="doc-list">
        <li><a href="#RandomVector-RandomVector">Constructor</a></li>
        <li>Properties</li>
        <ul class="doc-list">
            <li class="doc-list"><a href="#RandomVector-prop-type">type</a></li>
            <li class="doc-list"><a href="#RandomVector-prop-parameters">parameters</a></li>
            <li class="doc-list"><a href="#RandomVector-prop-dim">dim</a></li>
            <li class="doc-list"><a href="#RandomVector-prop-generator">generator</a></li>
        </ul>
        <li>Methods</li>
        <ul class="doc-list">
            <li class="doc-list"><a href="#RandomVector-method-concat">concat</a></li>
            <li class="doc-list"><a href="#RandomVector-method-vertcat">vertcat</a></li>
            <li class="doc-list"><a href="#RandomVector-method-getProbPolyhedron">getProbPolyhedron</a></li>
            <li class="doc-list"><a href="#RandomVector-method-cov">cov</a></li>
            <li class="doc-list"><a href="#RandomVector-method-mean">mean</a></li>
            <li class="doc-list"><a href="#RandomVector-method-pdf">pdf</a></li>
            <li class="doc-list"><a href="#RandomVector-method-getRealizations">getRealizations</a></li>
            <li class="doc-list"><a href="#RandomVector-method-plus">plus</a></li>
            <li class="doc-list"><a href="#RandomVector-method-mtimes">mtimes</a></li>
            <li class="doc-list"><a href="#RandomVector-method-horzcat">horzcat</a></li>
            <li class="doc-list"><a href="#RandomVector-method-uniform">uniform</a></li>
            <li class="doc-list"><a href="#RandomVector-method-gaussian">gaussian</a></li>
            <li class="doc-list"><a href="#RandomVector-method-exponential">exponential</a></li>
        </ul>
    </ul>
</ul>

{:#RandomVector}
### RandomVector
```
  Create a random vector object
  ==========================================================================
 
  Defines a random vector. We currently support:
  1. Gaussian    : Characterized by its mean vector and covariance matrix
  2. UserDefined : Characterized by a random number generator
 
 
  Usage:
  ------
 
  % Define a Gaussian random variable of mean 0 and standard deviation 2:
  GaussianRV = RandomVector('Gaussian', 0, 2^2);
 
  % Define a Gaussian random vector of mean [0;2] and covariance matrix 
  % eye(2):
  GaussianRV = RandomVector('Gaussian', [0;2], eye(2));
  % OR
  GaussianRV = RandomVector.gaussian([0;2], eye(2));
 
  % Define a beta-distributed 3-dimensional random vector with parameters
  % A=B=10
  BetaRV = RandomVector('UserDefined', @(N) betarnd(10,10,[3 N]));
    
  ==========================================================================
 
  RandomVector Properties:
  ------------------------
    type       - Random vector type (string)
    parameters - Random vector parameters (struct)
                    Stores mean and covariance for Gaussian random vector
    dim        - Random vector dimension (scalar)
    generator  - Random variable realization generator function (function 
                 handle); should take single numeric input and return an 
                 n x p matrix---where n is the dimension of the random vector
                 and p is the number of realizations (input value)
 
  RandomVector Methods:
  ---------------------
    RandomVector/RandomVector - Class constructor
    mean                      - Get the mean of random vector. If the mean is
                                not defined then an empirical mean is computed
                                using n_particle realizations [Default: 1e4]
    cov                       - Get the covariance of random vector. If the 
                                mean is not defined, then an empirical mean is 
                                computed using n_particle realizations. 
                                [Default: 1e4]
    pdf                       - Get the probability density function as an
                                anonymous function, defined only for Gaussian RV
    concat                    - Get the concatenated RV for a given time horizon
    getRealizations           - Generate realizations of the random vector
    getProbPolyhedron         - Get the probability of the random vector
                                lying in a user-specified polyhedron (MPT
                                object)
    Static
    ------
    gaussian                  - Get a Gaussian random vector for a
                                specified mean vector and covariance matrix
    exponential               - Get an exponential random vector for a
                                specified lambda vector
 
    Apart from these methods, for a RandomVector object rv, you can do:
    disp(rv)                  - Display information about rv
    F*rv, rv*F                - Both work as F x rv (NO TRANSPOSE enforced)
                                where F is an appropriately dimensioned matrix
    rv + v                    - Add v, a deterministic vector of appropriate
                                dimension, to rv
    [rv1;rv2;rv3]             - Concatenate multiple random vectors
  
  Notes:
  ------
  * MATLAB DEPENDENCY: Uses MATLAB's Statistics and Machine Learning Toolbox.
 
  
  =========================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
  

    Reference page in Doc Center
       doc RandomVector

```

{:#RandomVector-RandomVector}
### Constructor
```
  Create a random vector object
  ==========================================================================
 
  Defines a random vector. We currently support:
  1. Gaussian    : Characterized by its mean vector and covariance matrix
  2. UserDefined : Characterized by a random number generator
 
 
  Usage:
  ------
 
  % Define a Gaussian random variable of mean 0 and standard deviation 2:
  GaussianRV = RandomVector('Gaussian', 0, 2^2);
 
  % Define a Gaussian random vector of mean [0;2] and covariance matrix 
  % eye(2):
  GaussianRV = RandomVector('Gaussian', [0;2], eye(2));
  % OR
  GaussianRV = RandomVector.gaussian([0;2], eye(2));
 
  % Define a beta-distributed 3-dimensional random vector with parameters
  % A=B=10
  BetaRV = RandomVector('UserDefined', @(N) betarnd(10,10,[3 N]));
    
  ==========================================================================
 
  RandomVector Properties:
  ------------------------
    type       - Random vector type (string)
    parameters - Random vector parameters (struct)
                    Stores mean and covariance for Gaussian random vector
    dim        - Random vector dimension (scalar)
    generator  - Random variable realization generator function (function 
                 handle); should take single numeric input and return an 
                 n x p matrix---where n is the dimension of the random vector
                 and p is the number of realizations (input value)
 
  RandomVector Methods:
  ---------------------
    RandomVector/RandomVector - Class constructor
    mean                      - Get the mean of random vector. If the mean is
                                not defined then an empirical mean is computed
                                using n_particle realizations [Default: 1e4]
    cov                       - Get the covariance of random vector. If the 
                                mean is not defined, then an empirical mean is 
                                computed using n_particle realizations. 
                                [Default: 1e4]
    pdf                       - Get the probability density function as an
                                anonymous function, defined only for Gaussian RV
    concat                    - Get the concatenated RV for a given time horizon
    getRealizations           - Generate realizations of the random vector
    getProbPolyhedron         - Get the probability of the random vector
                                lying in a user-specified polyhedron (MPT
                                object)
    Static
    ------
    gaussian                  - Get a Gaussian random vector for a
                                specified mean vector and covariance matrix
    exponential               - Get an exponential random vector for a
                                specified lambda vector
 
    Apart from these methods, for a RandomVector object rv, you can do:
    disp(rv)                  - Display information about rv
    F*rv, rv*F                - Both work as F x rv (NO TRANSPOSE enforced)
                                where F is an appropriately dimensioned matrix
    rv + v                    - Add v, a deterministic vector of appropriate
                                dimension, to rv
    [rv1;rv2;rv3]             - Concatenate multiple random vectors
  
  Notes:
  ------
  * MATLAB DEPENDENCY: Uses MATLAB's Statistics and Machine Learning Toolbox.
 
  
  =========================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
  

    Reference page in Doc Center
       doc RandomVector

```

### Property: type
{:#RandomVector-prop-type}
```
  RandomVector/type
  ==================================================================
  
  String indicator of type of random vector
 
  Acceptable types:
    Gaussian
  
```

### Property: parameters
{:#RandomVector-prop-parameters}
```
  RandomVector/parameters
  ==================================================================
  
  Struct containing random vector parameter information; will vary for
  different random vector types. Currently only Gaussian random vector
  are supported.
 
  Gaussian type:
    parameters.mean       - Mean vector (p x 1)
    parameters.covariance - Covariance matrix (p x p)
  
```

### Property: dim
{:#RandomVector-prop-dim}
```
  RandomVector/dim
  ==================================================================
  
  Dimension of the random vector
  
```

### Property: generator
{:#RandomVector-prop-generator}
```
  RandomVector/generator
  ==================================================================
  
  Function to generate instances of the random variable
  
```

### Method: concat
{:#RandomVector-method-concat}
```
  Create a concatenated random vector of length time_horizon
  ====================================================================
  
  Inputs:
  -------
    obj           - RandomVector object (typically disturbance w_k)
    repeat_times  - The number of time steps of interest N
 
  Outputs:
  --------
    newobj        - RandomVector object (for disturbance, it can be used 
                    as W = [w_0^\top w_1^\top ... w_{N-1}^\top])
 
  Notes:
  ------
  * We make the independent and identical assumption to obtain the
    concatenated random vector
 
  ====================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```

### Method: vertcat
{:#RandomVector-method-vertcat}
```
  Vertical concatentation routine
  ====================================================================
  
  Inputs:
  -------
    obj             - A collection of RandomVector objects
 
  Outputs:
  --------
    newobj          - A RandomVector object that is the vertical
                      concatenation of random vectors
 
  Notes:
  ------
  * This function requires the objects that are being concatenated
    be of same RandomVector.type.
  * If the concatenation of a deterministic vector dv is desired with
    a RandomVector object rv, please use the following affine
    transformation-based command:
 
    new_rv = [dv; zeros(rv.dim,1)] + [zeros(length(dv), rv.dim);
                                      eye(rv.dim) ] * rv;
  
  ====================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
```

### Method: getProbPolyhedron
{:#RandomVector-method-getProbPolyhedron}
```
  Compute the probability of a random vector lying in a polyhedron 
  ====================================================================
  Probability computation is done via Monte-Carlo (quasi Monte-Carlo in
  Gaussian case). A distribution-independent lower bound on the number
  of particles needed is obtained via Hoeffding's inequality.
 
  Hoeffding's inequality states that 
 
        Prob{|X-E[X]|\geq \delta} \leq \beta, with
         \beta = 2 * exp(-2 * n_particles * \delta^2). 
 
  Here, X is a Bernoulli random variable, \delta is the
  desired_accuracy, \beta is the failure risk (the probability of the
  statement |X-E[X]|\geq \delta fails). 
 
  In this case, X corresponds to the indicator function of the polytope
  composed with the random vector. Given a failure risk \beta and
  desired accuracy \delta, we backcompute the number of particles,
  
        n_particles = -ln(\beta/2) / (2 * delta^2) 
  
  enforces the condition that empirical mean does not deviate more than
  delta from the true mean, by Hoeffding's inequality.
  
  Inputs:
  -------
    obj             - RandomVector object
    test_polyhedron - Polyhedron object (polytope whose probability of
                      occurrence is of interest)
    desired_accuracy- [Optional] Maximum absolute deviation from the
                      true probability estimate [Default: 1e-2]
 
  Outputs:
  --------
    covar           - Probability of the random vector lying in the
                      given polytope
  Notes:
  ------
  * Due to the inverse-square dependence on the desired_accuracy, we
    impose a hard lower bound of 1e-2 on the desired_accuracy. This
    leads to the requirement of 2e5 particles.
  * We set the default failure risk as 2e-15.
  * Ill-formed (mean/cov has Inf/NaN) Gaussian random vectors return
    zero probability.
 
  ====================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```

### Method: cov
{:#RandomVector-method-cov}
```
  Convenience method for accessing the covariance of a random vector
  ====================================================================
  
  Inputs:
  -------
    obj         - RandomVector object
    n_particles - [Optional] Number of particles to use in sample 
                  covariance estimation [Default: 1e6]
 
  Outputs:
  --------
    covar       - Covariance of random vector
 
  Notes:
  ------
  * For Random Vectors of type 'Gaussian', we return the exact
    covariance matrix
  * For Random Vectors of type 'UserDefined', we return MATLAB's
    estimated covariance matrix obtained from n_particles samples
  
  ====================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```

### Method: mean
{:#RandomVector-method-mean}
```
  Convenience method for accessing the (sample) mean of a random vector
  ====================================================================
  
  Inputs:
  -------
    obj         - RandomVector object
    n_particles - [Optional] Number of particles to use in sample 
                  mean estimation [Default: 1e6]
 
 
  Outputs:
  --------
    m           - Mean of random vector
 
  Notes:
  ------
  * For Random Vectors of type 'Gaussian', we return the exact mean
  * For Random Vectors of type 'UserDefined', we return MATLAB's
    estimated mean obtained from n_particles samples
  
  ====================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```

### Method: pdf
{:#RandomVector-method-pdf}
```
  Get the pdf of a random vector
  ====================================================================
  
  Inputs:
  -------
    obj - RandomVector object
 
  Outputs:
  --------
    pdf - Probability density function (anonymous function handle)
 
  Notes:
  ------
  * This code is not tested
  * Only Random Vectors of type 'Gaussian' currently have a pdf.
  * Other random vector types will throw an error
  * The anonymous function used for the definition of obj.pdf transposes 
    the accepted column vector for using mvnpdf.
  * RandomVector.pdf takes in arguments of the form N_points x
    random_vector_dim
  
  ====================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```

### Method: getRealizations
{:#RandomVector-method-getRealizations}
```
  Generate n_realizations realizations of the random vector
  ====================================================================
  
  Inputs:
  -------
    obj            - RandomVector object (typically disturbance w_k)
    n_realizations - Number of realizations desired
 
  Outputs:
  --------
    xs             - Realizations matrix of dimension 
                     obj.dim x n_realizations. Each realization is given 
                     as a column vector.
 
  Notes:
  ------
  * In case of Gaussian random vectors, mvnrnd is used
  * In case of UserDefined random vectors, the user-provided
    generator is used
 
  ====================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```

### Method: plus
{:#RandomVector-method-plus}
```
  Override of MATLAB plus command
  ====================================================================
  
  Inputs:
  -------
    obj - RandomVector object
    v   - Deterministic vector to be added to the random vector OR
          a RandomVector object
 
  Outputs:
  --------
    newobj - RandomVector object (obj + v)
 
  Notes:
  ------
  * While this function updates the generator for a UserDefined random
    vector, it is HIGHLY RECOMMENDED to redefine the random vector
    separately with an updated generator function to avoid nested
    generator functions
  ====================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```

### Method: mtimes
{:#RandomVector-method-mtimes}
```
  Override of MATLAB multiplication command
  ====================================================================
  
  Inputs:
  -------
    obj - RandomVector object
    F   - Linear transformation matrix for multiplication
 
  Outputs:
  --------
    newobj - RandomVector object (F*obj)
 
  Notes:
  ------
  * While this function updates the generator for a UserDefined random
    vector, it is HIGHLY RECOMMENDED to redefine the random vector
    separately with an updated generator function to avoid nested
    generator functions
  ====================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```

### Method: horzcat
{:#RandomVector-method-horzcat}
```
  Horizontal concatenation prevention routine!
  ====================================================================
  
  Notes:
  ------
  * This function just throws an error since horizontal
    concatenation produces random matrix!
  
  ====================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```

### Method: uniform
{:#RandomVector-method-uniform}
```
RandomVector.uniform is a function.
    rv = uniform(lb, ub)
```

### Method: gaussian
{:#RandomVector-method-gaussian}
```
  Convenience method for creating Gaussian random vectors
  ====================================================================
  
  Inputs:
  -------
    mu     - Gaussian mean vector (column vector)
    covar  - Gaussian covariance matrix (square matrix)
 
  Outputs:
  --------
    rv     - RandomVector object
 
  ====================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```

### Method: exponential
{:#RandomVector-method-exponential}
```
  Convenience method for creating exponential random vectors
  ====================================================================
  
  Static method used to conveniently create n-dimensional exponential
  random vectors. Vectors are generated using MATLAB's exprnd function.
  
  ====================================================================
  
  Inputs:
  -------
    mu  - Exponential mean | Must be a column vector
 
  Outputs:
  --------
    rv - RandomVector object
 
  ====================================================================
  
  This function is part of the Stochastic Reachability Toolbox.
  License for the use of this function is given in
       https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  
 
```

