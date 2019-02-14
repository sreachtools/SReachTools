---
layout: docs
title: qscmvnv.m
---

```
   [ P E ] = QSCMVNV( M, R, A, CN, B )
     uses a randomized quasi-random rule with m points to estimate an
     MVN probability for positive semi-definite covariance matrix r,
     with constraints a < cn*x < b. If r is nxn and cn is kxn, then
     a and b must be column k-vectors.
    Probability p is output with error estimate e.
     Example use:
      r = [ 4 3 2 1; 3 5 -1 1; 2 -1 4 2; 1 1 2 5 ];
      a = [ -inf 1 -5 ]'; b = [ 3 inf 4 ]';
      cn = [ 1 2 3 -2; 2 4 1 2; -2 3 4 1 ];
      [ p e ] = qscmvnv( 5000, r, a, cn, b ); disp([ p e ])
 
   This function uses an algorithm given in the paper by Alan Genz:
    "Numerical Computation of Multivariate Normal Probabilities", in
      J. of Computational and Graphical Stat., 1(1992), 141-149.
   The primary references for the numerical integration are 
    "On a Number-Theoretical Integration Method"
      H. Niederreiter, Aequationes Mathematicae, 8(1972), 304-11, and
    "Randomization of Number Theoretic Methods for Multiple Integration"
      R. Cranley and T.N.L. Patterson, SIAM J Numer Anal, 13(1976), 904-14.
 
    Alan Genz is the author of this function and following Matlab functions.
           Alan Genz, WSU Math, PO Box 643113, Pullman, WA 99164-3113
           Email : AlanGenz@wsu.edu
 
 
  Copyright (C) 2014, Alan Genz,  All rights reserved.               
 
  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided the following conditions are met:
    1. Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer.
    2. Redistributions in binary form must reproduce the above copyright
       notice, this list of conditions and the following disclaimer in 
       the documentation and/or other materials provided with the 
       distribution.
    3. The contributor name(s) may not be used to endorse or promote 
       products derived from this software without specific prior 
       written permission.
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS 
  OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR 
  TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```
