function [A_norm, b_norm] = normalizeForParticleControl(A,b)
% Normalize (A,b) to (A_norm, b_norm) where b_norm \in ((\infty,1e-3) U {1})^L
% ============================================================================
% 
% Given a set {z: Az<= b}, this function computes (A_norm, b_norm) such that the
% same set is represented by {z: A_norm z<= b_norm} and b_norm \in
% ((\infty,1e-3) U {1})^L. The advantage with this normalization is that the
% bigM used in particle control can be significantly smaller (100 is sufficient)
% in most cases, since the value taken by b is either <= 0 or 1. 
%
% See also SReachPointVoA.m
%
% =============================================================================
% [A_norm, b_norm] = normalizeForParticleControl(A,b)
%
% Inputs:
% -------
%   (A,b) - Half space representation of a given polytope
%
% Outputs:
% --------
%   (A_norm, b_norm)  
%         - Normalized Half space representation of a given polytope
%
% ============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://sreachtools.github.io/license/
% 
%

    zero_defn = 1e-3;
    A_norm = A;
    b_norm = b;
    
    % Normalize b that are above zero_defn
    % disp(nnz(b >= zero_defn));
    b_norm(b >= zero_defn) = 1;
    A_norm(b >= zero_defn, :) = A(b >= zero_defn, :)./b(b >= zero_defn);    
end
