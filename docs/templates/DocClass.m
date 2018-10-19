classdef DocClass
%  One-sentence description of the class
% ===========================================================================
%
% More detailed description of the class. This area can have a minor paragraph
% structure but should still be succint description of the class. This class
% is an empty placeholder class used to demonstrate class documentation style.
%
% Usage:
% ------
% % Self-contained examples
%
% sys = DocClass('blue', 'suede', 'shoes', 'optionName', 4);
%   
% ===========================================================================
%
% DocClass Properties:
% ---------------------
%   property_1 - First class property (type if applicable)
%   property_2 - Second class property (type if applicable)
%   property_3 - Third class property (type if applicable)
% 
% See also <<any applicable reference functions or classes>>
% 
% DocClass Methods:
% ------------------
%   DocClass/DocClass - Constructor
%   someRandomMethod  - A dummy method that returns a random number between 
%                       zero and one
%
% ===========================================================================
%
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
%
%
    properties
        property_1

        % some more detailing about property 2 because it's the best
        property_2
        property_3
    end

    methods
        function obj = DocClass(in1, in2, in3, varargin)
        %  DocClass constructor
        % ====================================================================
        %
        % Example class constructor. Follow documentation style guidelines for 
        % standard function. See doc/templates/standardFunction.m
        %
        % Usage:
        % ------
        % sys = DocClass('blue', 'suede', 'shoes', 'optionName', 4);
        %
        % ====================================================================
        %
        % obj = DocClass(in1, in2, in3, Name Value)
        % 
        % Inputs:
        % -------
        %   in1 - Input 1
        %   in2 - Input 2
        %   in3 - Input 3
        % 
        %   -----------------------------------------------------
        %   Name              | Value
        %   -----------------------------------------------------
        %   Anything          | Anything esle
        % 
        % See also Polyhedron, StochasticDisturbance
        % 
        % Outputs:
        % --------
        % obj - DocClass object
        %
        % Notes:
        % ------
        % * Put any useful notes here
        % 
        % ====================================================================
        % 
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        % 
            obj.property_1 = in1;
            obj.property_2 = in2;
            obj.property_3 = in3;
        end

        function a = someRandomMethod(obj)
        %  Random method to return a random
        % number
        % ====================================================================
        %
        % Dummy method demonstrating an example of a class method. Returns a 
        % random number in [0,1]. Follow documentation style guidelines for 
        % standard function. See doc/templates/standardFunction.m
        %
        % Usage:
        % ------
        % sys = DocClass('blue', 'suede', 'shoes', 'optionName', 4);
        % disp(sys.someRandomMethod())
        %
        % ====================================================================
        %
        % a = someRandomMethod()
        % 
        % Inputs: None (Do not list obj)
        % 
        % Outputs:
        % --------
        % a - Random number (double)
        %
        % ====================================================================
        % 
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        % 
            a = rand();
        end
    end
end