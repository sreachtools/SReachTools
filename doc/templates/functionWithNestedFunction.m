function [a, b, c] = functionWithNestedFunction(d, e, f)
%  Function with nested function 
% documentation example
% ============================================================================
%
% Example of standard function documentation
%
% Usage:
% ------
% % so easy
% [a, b, c] = functionWithNestedFunction(1, 2, 3)
%
% ============================================================================
%
% [a, b, c] = functionWithNestedFunction(d, e, f)
% 
% Inputs:
% -------
%   d - Input
%   e - Input
%   f - Input
% 
% See also <<other useful functions relating to this>>
% 
% Outputs:
% --------
%   a - Output
%   b - Output
%   c - Output
%
% Notes:
% ------
% * a = 1
% * b = 2
% * c = 3
% 
% ============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 
% 

    a = 1;
    b = 2;
    c = 3;
    
end

function output = nestedFunction(in1, in2)
%  Nested function documentation example
% ============================================================================
%
% Nested function documentation is very similar to standard function
% documentation, however since it will only be seen by a developer some details
% are not necessary. Specifically, the usage is not needed since the function
% cannot be explicitly called from the MATLAB console. Simply list its
% appropriate parent function.
%
% Usage:
% ------
% Nested function of functionWithNestedFunction
%
% ============================================================================
%
% output = nestedFunction(in1, in2)
% 
% Inputs:
% -------
%   in1 - Input
%   in2 - Input
% 
% See also <<other useful functions relating to this>>
% 
% Outputs:
% --------
% output - Output
% 
% ============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 
% 

    output = 'super-awesome!';
end