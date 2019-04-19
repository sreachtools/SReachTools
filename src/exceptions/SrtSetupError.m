classdef SrtSetupError < SrtBaseException
% Custom exception object for SReachTools setup errors
% ============================================================================
% 
% Customized class for generating SReachTools setup errors, subclass of the 
% standard MATLAB SrtBaseException class
%
% Usage:
% ------
% exc = SrtSetupError('error message')
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
        mnemonic = 'setupError';
    end
    
    methods
        function obj = SrtSetupError(varargin)
            obj@SrtBaseException(SrtSetupError.mnemonic, varargin{:});
       end
    end
    
    methods (Static)
        function id = getErrorId()
            id = [SrtBaseException.getErrorComponent(), ':', ...
                SrtInvalidArgsError.mnemonic];
        end
    end
end

