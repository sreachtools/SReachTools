% Stochastic Reachability Toolbox (SReachTools)
% Version 1.2.27 (R2017b) 5-December-2018
% 
% The Stochastic Reachability Toolbox (SReachTools) is an open-source MATLAB 
% toolbox for performing stochastic verification and reachability analysis.
%
% For questions and comments, visit the project's website at
%
%    https://unm-hscl.github.io/SReachTools/
%
% Compute the set of safe initial states and associated admissible controllers 
% to drive the system to stay with in a target tube (a time-varying collection 
% of target sets). The toolbox can handle linear (time-invariant or
% time-varying) systems with additive Gaussian.
%     - Stochastic reach set (approximations) via 
%           1. Lagrangian methods
%           2. Convex chance constrained optimization
%           3. Fourier transform (Genz's algorithm and MATLAB's patternsearch)
%     - Optimal controller synthesis
%           1. Open-loop controller synthesis
%               1. Convex chance constrained optimization
%               2. Fourier transform (Genz's algorithm and MATLAB's
%                  patternsearch)
%               3. Particle control
%               4. Voronoi-based undersampled particle control
%           2. Affine controller synthesis
%               1. Difference-of-convex chance constrained optimization
%               2. Voronoi-based undersampled particle control
%     - Dynamic programming methods for 'exact' computation of the reach
%       probability using backward recursion
%     - Forward stochastic reachability analysis via Gaussian random vector
%       properties
%     
% Please see examples/ folder for demonstration of various features of
% SReachTools.
%
% Documentation of this toolbox is available online at
% https://unm-hscl.github.io/SReachTools/docs/. Further, use <help fun_name>
% to get the MATLAB docstring associated with a function fun_name.

%   docs2md - Script to update documentation
%   srtinit - Initialization function



