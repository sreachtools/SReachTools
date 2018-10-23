% Stochastic Reachability Toolbox (SReachTools)
% Version 1.0.0 (R2017b) 19-October-2018
% 
% The Stochastic Reachability Toolbox (SReachTools) is an open-source MATLAB 
% toolbox for performing stochastic verification and reachability analysis.
%
% For questions and comments, visit the project's website at
%
%    https://unm-hscl.github.io/SReachTools/
%
% Basic Functionality: Compute the set of safe initial states and associated
%   admissible controllers to drive the system to stay with in a target tube (a
%   time-varying collection of target tubes)
%     - Stochastic reach set (approximations) via 
%           1. Lagrangian methods
%           2. Convex chance constrained optimization
%           3. Fourier transform (Genz's algorithm and MATLAB's patternsearch)
%     - Optimal controller synthesis
%           1. Open-loop controller synthesis
%               1. Convex chance constrained optimization
%               2. Fourier transform (Genz's algorithm and MATLAB's
%                  patternsearch)
%               3. Particle filter
%           1. Affine controller synthesis
%               1. Difference of convex and chance constrained optimization
%     - Dynamic programming methods for 'exact' computation of the reach
%       probability using backward recursion
%     - Forward stochastic reachability analysis via Gaussian random vector
%       properties
%     

