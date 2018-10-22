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
            <li class="doc-list"><a href="#RandomVector-prop-pdf">pdf</a></li>
        </ul>
        <li>Methods</li>
        <ul class="doc-list">
        </ul>
    </ul>
</ul>

{:#RandomVector}
### RandomVector
```
  Create a random vector object
  ==========================================================================
 
  Defines a random vector with a standard probability density function (pdf)
 
  We support the following pdfs:
      1. Gaussian
 
  Usage:
  ------
 
  % Define a Gaussian random variable of mean 0 and standard deviation 2:
  GaussianRV = RandomVector('Gaussian', ...
                            0, ...
                            2^2);
 
  % Define a Gaussian random vector of mean [0;2] and covariance matrix 
  % eye(2):
  GaussianRV = RandomVector('Gaussian', ...
                            [0;2], ...
                            eye(2));
    
  ==========================================================================
 
  RandomVector Properties:
  ------------------------
    type       - Random vector type (string)
    parameters - System parameters (struct)
    dim        - Random vector dimension (scalar)
    pdf        - Probability density function (function handle)
 
  RandomVector Methods:
  ---------------------
    RandomVector/RandomVector - Class constructor
  
  Notes:
  ------
  * MATLAB DEPENDENCY: Uses MATLAB's Statistics and Machine Learning Toolbox.
                       Needs mvnpdf
  * Currently only supports Gaussian random vectors
  * Requires the mean and the covariance matrices to be non-empty column
    vector and a symmetric matrix respectively
  * The anonymous function used for the definition of obj.pdf transposes the
    accepted column vector for using mvnpdf.
  * RandomVector.pdf takes in arguments of the form N_points x
    random_vector_dim
  
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
 
  Defines a random vector with a standard probability density function (pdf)
 
  We support the following pdfs:
      1. Gaussian
 
  Usage:
  ------
 
  % Define a Gaussian random variable of mean 0 and standard deviation 2:
  GaussianRV = RandomVector('Gaussian', ...
                            0, ...
                            2^2);
 
  % Define a Gaussian random vector of mean [0;2] and covariance matrix 
  % eye(2):
  GaussianRV = RandomVector('Gaussian', ...
                            [0;2], ...
                            eye(2));
    
  ==========================================================================
 
  RandomVector Properties:
  ------------------------
    type       - Random vector type (string)
    parameters - System parameters (struct)
    dim        - Random vector dimension (scalar)
    pdf        - Probability density function (function handle)
 
  RandomVector Methods:
  ---------------------
    RandomVector/RandomVector - Class constructor
  
  Notes:
  ------
  * MATLAB DEPENDENCY: Uses MATLAB's Statistics and Machine Learning Toolbox.
                       Needs mvnpdf
  * Currently only supports Gaussian random vectors
  * Requires the mean and the covariance matrices to be non-empty column
    vector and a symmetric matrix respectively
  * The anonymous function used for the definition of obj.pdf transposes the
    accepted column vector for using mvnpdf.
  * RandomVector.pdf takes in arguments of the form N_points x
    random_vector_dim
  
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

### Property: pdf
{:#RandomVector-prop-pdf}
```
  RandomVector/pdf
  ==================================================================
  
  Probability density function of random vector
  
```

