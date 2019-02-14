function [ p, e ] = qscmvnv( m, r, a, cn, b )
% 
%  [ P E ] = QSCMVNV( M, R, A, CN, B )
%    uses a randomized quasi-random rule with m points to estimate an
%    MVN probability for positive semi-definite covariance matrix r,
%    with constraints a < cn*x < b. If r is nxn and cn is kxn, then
%    a and b must be column k-vectors.
%   Probability p is output with error estimate e.
%    Example use:
%     r = [ 4 3 2 1; 3 5 -1 1; 2 -1 4 2; 1 1 2 5 ];
%     a = [ -inf 1 -5 ]'; b = [ 3 inf 4 ]';
%     cn = [ 1 2 3 -2; 2 4 1 2; -2 3 4 1 ];
%     [ p e ] = qscmvnv( 5000, r, a, cn, b ); disp([ p e ])
%
%  This function uses an algorithm given in the paper by Alan Genz:
%   "Numerical Computation of Multivariate Normal Probabilities", in
%     J. of Computational and Graphical Stat., 1(1992), 141-149.
%  The primary references for the numerical integration are 
%   "On a Number-Theoretical Integration Method"
%     H. Niederreiter, Aequationes Mathematicae, 8(1972), 304-11, and
%   "Randomization of Number Theoretic Methods for Multiple Integration"
%     R. Cranley and T.N.L. Patterson, SIAM J Numer Anal, 13(1976), 904-14.
%
%   Alan Genz is the author of this function and following Matlab functions.
%          Alan Genz, WSU Math, PO Box 643113, Pullman, WA 99164-3113
%          Email : AlanGenz@wsu.edu
%
%
% Copyright (C) 2014, Alan Genz,  All rights reserved.               
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided the following conditions are met:
%   1. Redistributions of source code must retain the above copyright
%      notice, this list of conditions and the following disclaimer.
%   2. Redistributions in binary form must reproduce the above copyright
%      notice, this list of conditions and the following disclaimer in 
%      the documentation and/or other materials provided with the 
%      distribution.
%   3. The contributor name(s) may not be used to endorse or promote 
%      products derived from this software without specific prior 
%      written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
% FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
% COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
% BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS 
% OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR 
% TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%
% Initialization
%
[ as ch bs clg n ] = chlsrt( r, a, cn, b );
ci = Phi(as(1)); dci = Phi(bs(1)) - ci; p = 0; e = 0; 
ns = 12; nv = fix( max( [ m/ns 1 ] ) ); 
ps = sqrt(primes(5*(n+1)*log(n+1)/4)); q = ps(1:n-1)'; % Richtmyer generators
%
% Randomization loop for ns samples
%
for i = 1 : ns
  % periodizing transformation 
  xx(:,1:nv) = abs( 2*mod( q*[1:nv] + rand(n-1,1)*ones(1,nv), 1 ) - 1 );
  vp =   mvndnv( n, as, ch, bs, clg, ci, dci, xx, nv ); 
   d = ( mean(vp) - p )/i; p = p + d;
  if abs(d) > 0, e = abs(d)*sqrt( 1 + ( e/d )^2*(i-2)/i );
  else, if i > 1, e = e*sqrt( ( i - 2 )/i ); end
  end
end, e = 3*e; % error estimate is 3 x standard error with ns samples.
%
% end qscmvnv
%
function p = mvndnv( n, a, ch, b, clg, ci, dci, x, nv )
%
%  Transformed integrand for computation of MVN probabilities. 
%
y = zeros(n-1,nv); on = ones(1,nv); 
c = ci*on; dc = dci*on; p = dc; li = 2; lf = 1;
for i = 2 : n, y(i-1,:) = Phinv( c + x(i-1,:).*dc ); lf = lf + clg(i); 
   if lf < li, c = 0; dc = 1;
   else, s = ch(li:lf,1:i-1)*y(1:i-1,:); 
     ai =          max( max( a(li:lf)*on - s, [], 1 ), -36 ); 
     bi = max( ai, min( min( b(li:lf)*on - s, [], 1 ),   9 ) ); 
     c = Phi(ai); dc = Phi(bi) - c; p = p.*dc; 
   end, li = li + clg(i);
end 
%
% end mvndnv
%
function [ ap, ch, bp, clg, np ] = chlsrt( r, a, cn, b )
%
%  Computes permuted lower Cholesky factor ch for covariance r which 
%   may be singular, combined with contraints a < cn*x < b, to
%   form revised lower triangular constraint set ap < ch*x < bp; 
%   clg contains information about structure of ch: clg(1) rows for 
%   ch with 1 nonzero, ..., clg(np) rows with np nonzeros.
%
ep = 1e-10; % singularity tolerance;
%
[n n] = size(r); y = zeros(n,1); clg = zeros(1,n); ap = a; bp = b; 
[m n] = size(cn); ch = cn; c = r; d = sqrt(max(diag(c),0));
for i = 1 : n, di = d(i);
  if di>0, c(:,i) = c(:,i)/di; c(i,:) = c(i,:)/di; ch(:,i) = ch(:,i)*di; end
end
%
%  determine r factors and form revised constraint matrix ch
%
[V,D] = eig(c); [d,pm] = sort(diag(D)','descend'); d = sqrt(max(d,0));  
np = sum(d>0); c = V(:,pm(1:np)).*( ones(n,1)*d(1:np) ); ch = ch*c; 
%
% use right reflectors to reduce ch to lower triangular
%
for i = 1 : min( np-1, m ), epi = ep*i; vm = 1; lm = i;
  %
  % permute rows so that smallest variance variables are first.
  %
  for l = i : m, v = ch(l,1:np); s = v(1:i-1)*y(1:i-1); 
    ss = max( sqrt( sum( v(i:np).^2 ) ), epi ); 
    al = ( ap(l) - s )/ss; bl = ( bp(l) - s )/ss; 
    dna = 0; dsa = 0; dnb = 0; dsb = 1;
    if al > -9, dna = phi(al); dsa = Phi(al); end
    if bl <  9, dnb = phi(bl); dsb = Phi(bl); end
    if dsb - dsa > epi
      if      al <= -9, mn =      -dnb; vr =         -bl*dnb;
      elseif  bl >=  9, mn = dna;       vr = al*dna; 
      else,             mn = dna - dnb; vr = al*dna - bl*dnb; 
      end, mn = mn/( dsb - dsa ); vr = 1 + vr/( dsb - dsa ) - mn^2;
    else, mn = ( al + bl )/2; vr = 0;
      if al <= -9, mn = bl; elseif bl >=  9, mn = al; end, 
    end, if vr <= vm, lm = l; vm = vr; y(i) = mn; end
  end, v = ch(lm,1:np);
  if lm > i, ch([i lm],1:np) = ch([lm i],1:np); 
    ap([i lm]) = ap([lm i]); bp([i lm]) = bp([lm i]); 
  end, ch(i,i+1:np) = 0; ss = sum( v(i+1:np).^2 );
  if ss > epi, ss = sqrt( ss + v(i)^2 ); if v(i) < 0, ss = -ss; end
    ch(i,i) = -ss; v(i) = v(i) + ss; vt = v(i:np)'/( ss*v(i) );
    ch(i+1:m,i:np) = ch(i+1:m,i:np) - ch(i+1:m,i:np)*vt*v(i:np); 
  end
end
%
% scale and sort constraints
%
for i = 1 : m, v = ch(i,1:np); clm(i) = min(i,np); 
  jm = 1; for j = 1 : clm(i), if abs(v(j)) > ep*j, jm = j; end, end 
  if jm < np, v(jm+1:np) = 0; end, clg(jm) = clg(jm) + 1; 
  at = ap(i); bt = bp(i); j = i;
  for l = i-1 : -1 : 1, if jm >= clm(l), break, end
    ch(l+1,1:np) = ch(l,1:np); j = l;
    ap(l+1) = ap(l); bp(l+1) = bp(l); clm(l+1) = clm(l);
  end
  clm(j) = jm; vjm = v(jm); ch(j,1:np) = v/vjm; 
  ap(j) = at/vjm; bp(j) = bt/vjm;
  if vjm < 0, tl = ap(j); ap(j) = bp(j); bp(j) = tl; end
end
j = 0; for i = 1 : np, if clg(i) > 0, j = i; end, end, np = j;
%
% combine constraints for first variable
%
if clg(1) > 1 
  ap(1) = max( ap(1:clg(1)) ); bp(1) = max( ap(1), min( bp(1:clg(1)) ) ); 
  ap(2:m-clg(1)+1) = ap(clg(1)+1:m); bp(2:m-clg(1)+1) = bp(clg(1)+1:m);
  ch(2:m-clg(1)+1,:) = ch(clg(1)+1:m,:); clg(1) = 1;
end
%
% end chlsrt
%
%  Standard statistical normal distribution pdf, cdf, and inverse
function p = phi(z),   p = normpdf(z);         % p = exp(-z.^2/2)/sqrt(2*pi);
function p = Phi(z),   p = normcdf(z);         % p = erfc( -z/sqrt(2) )/2;  
function z = Phinv(w), z = norminv(w);         % z = -sqrt(2)*erfcinv( 2*w );

