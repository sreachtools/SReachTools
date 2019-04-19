classdef SrtInternalError < SrtBaseException
% Custom exception object for SReachTools internal errors
% ============================================================================
% 
% Customized class for generating SReachTools internal errors, subclass of the 
% standard MATLAB SrtBaseException class
%
% Usage:
% ------
% exc = SrtInternalError('error message')
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
        mnemonic = 'internal';
    end
    
    methods
        function obj = SrtInternalError(varargin)
            obj@SrtBaseException(SrtInternalError.mnemonic, varargin{:}); 
        end
    end

    methods (Static)
        function id = getErrorId()
            id = [SrtBaseException.getErrorComponent(), ':', ...
                SrtInternalError.mnemonic];
        end
    end
    
end

