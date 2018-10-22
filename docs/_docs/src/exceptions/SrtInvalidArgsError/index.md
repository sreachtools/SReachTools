---
layout: docs
title: SrtInvalidArgsError.m
---

<ul class="doc-list">
    <li class="doc-list"><a href="#SrtInvalidArgsError">SrtInvalidArgsError</a></li>
    <ul class="doc-list">
        <li><a href="#SrtInvalidArgsError-SrtInvalidArgsError">Constructor</a></li>
        <li>Properties</li>
        <ul class="doc-list">
        </ul>
        <li>Methods</li>
        <ul class="doc-list">
            <li class="doc-list"><a href="#SrtInvalidArgsError-method-withFunctionName">withFunctionName</a></li>
            <li class="doc-list"><a href="#SrtInvalidArgsError-method-getErrorId">getErrorId</a></li>
        </ul>
    </ul>
</ul>

{:#SrtInvalidArgsError}
### SrtInvalidArgsError
```
  Custom exception object for SReachTools invalid arguments errors
  ============================================================================
  
  Customized class for generating SReachTools invalid arguments errors, subclass 
  of the standard SrtBaseException class
 
  Usage:
  ------
  exc = SrtInvalidArgsError('error message')
 
  ============================================================================
 
  See also MException
 
  ============================================================================
 
    This function is part of the Stochastic Optimal Control Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  

    Reference page in Doc Center
       doc SrtInvalidArgsError

```

{:#SrtInvalidArgsError-SrtInvalidArgsError}
### Constructor
```
  Custom exception object for SReachTools invalid arguments errors
  ============================================================================
  
  Customized class for generating SReachTools invalid arguments errors, subclass 
  of the standard SrtBaseException class
 
  Usage:
  ------
  exc = SrtInvalidArgsError('error message')
 
  ============================================================================
 
  See also MException
 
  ============================================================================
 
    This function is part of the Stochastic Optimal Control Toolbox.
    License for the use of this function is given in
         https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
  

    Reference page in Doc Center
       doc SrtInvalidArgsError

```

### Method: withFunctionName
{:#SrtInvalidArgsError-method-withFunctionName}
```
  Throw invalid args and provide function name that received the 
  invalid call
  =====================================================================
  
  Method to preformat error strings to throw to the user the specific
  function that was used when passing the invalid args error. Will 
  display message:
  
        Invalid arguments provided to << function name >>
  
  =====================================================================
  
```

### Method: getErrorId
{:#SrtInvalidArgsError-method-getErrorId}
```
SrtInvalidArgsError.getErrorId is a function.
    id = getErrorId()
```

