---
layout: docs
title: allcomb.m
---

```
  ALLCOMB - All combinations
     B = ALLCOMB(A1,A2,A3,...,AN) returns all combinations of the elements
     in the arrays A1, A2, ..., and AN. B is P-by-N matrix where P is the product
     of the number of elements of the N inputs. 
     This functionality is also known as the Cartesian Product. The
     arguments can be numerical and/or characters, or they can be cell arrays.
 
     Examples:
        allcomb([1 3 5],[-3 8],[0 1]) % numerical input:
        % -> [ 1  -3   0
        %      1  -3   1
        %      1   8   0
        %        ...
        %      5  -3   1
        %      5   8   1 ] ; % a 12-by-3 array
 
        allcomb('abc','XY') % character arrays
        % -> [ aX ; aY ; bX ; bY ; cX ; cY] % a 6-by-2 character array
 
        allcomb('xy',[65 66]) % a combination -> character output
        % -> ['xA' ; 'xB' ; 'yA' ; 'yB'] % a 4-by-2 character array
 
        allcomb({'hello','Bye'},{'Joe', 10:12},{99999 []}) % all cell arrays
        % -> {  'hello'  'Joe'        [99999]
        %       'hello'  'Joe'             []
        %       'hello'  [1x3 double] [99999]
        %       'hello'  [1x3 double]      []
        %       'Bye'    'Joe'        [99999]
        %       'Bye'    'Joe'             []
        %       'Bye'    [1x3 double] [99999]
        %       'Bye'    [1x3 double]      [] } ; % a 8-by-3 cell array
 
     ALLCOMB(..., 'matlab') causes the first column to change fastest which
     is consistent with matlab indexing. Example: 
       allcomb(1:2,3:4,5:6,'matlab') 
       % -> [ 1 3 5 ; 1 4 5 ; 1 3 6 ; ... ; 2 4 6 ]
 
     If one of the N arguments is empty, ALLCOMB returns a 0-by-N empty array.
     
     See also NCHOOSEK, PERMS, NDGRID
          and NCHOOSE, COMBN, KTHCOMBN (Matlab Central FEX)
```



## License
Available at `./src/helperFunctions/allcomb_license.txt`.

Copyright (c) 2018, Jos (10584)
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in
the documentation and/or other materials provided with the distribution

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.