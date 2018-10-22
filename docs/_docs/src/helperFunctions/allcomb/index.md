---
layout: docs
title: allcomb.m
---

```
  ALLCOMB - All combinations
     B = ALLCOMB(A1,A2,A3, ...,AN) returns all combinations of the elements
     in the arrays A1, A2, ..., and AN. B is P-by-N matrix is which P is the product
     of the number of elements of the N inputs. This functionality is also
     known as the Cartesian Product. The arguments can be numerical and/or
     characters, or they can be cell arrays.
 
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
 
        allcomb('xy',[65 66]) % a combination
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
 
     If one of the arguments is empty, ALLCOMB returns a 0-by-N empty array.
     
     See also NCHOOSEK, PERMS, NDGRID
          and NCHOOSE, COMBN, KTHCOMBN (Matlab Central FEX)
```
