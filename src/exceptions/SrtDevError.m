classdef SrtDevError < SrtBaseException
% Custom exception object for SReachTools development errors
% ============================================================================
% 
% Customized class for generating SReachTools development errors, subclass of the 
% standard MATLAB SrtBaseException class
%
% Usage:
% ------
% exc = SrtDevError('error message')
%
% ============================================================================
%
% See also MException
%
% ============================================================================
%
%   This function is part of the Stochastic Optimal Control Toolbox.
%   License for the use of this function is given in
%        https://sreachtools.github.io/license/
% 
    
    properties (Constant, Access = private)
        mnemonic = 'runtime';
    end
    
    methods
        function obj = SrtDevError(varargin)
            obj@SrtBaseException(SrtDevError.mnemonic, varargin{:}); 
            obj.message = sprintf([obj.message, ...
                 '\nThis might be a bug in the SReachTools. Please ', ...
                 'contact the devs at ', ...
                 'https://groups.google.com/d/forum/sreachtools with a ', ...
                 'minimal working example. Thank you!']);
        end
    end

    methods (Static)
        function id = getErrorId()
            id = [SrtBaseException.getErrorComponent(), ':', ...
                SrtDevError.mnemonic];
        end
    end
end

