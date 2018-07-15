classdef SrtSetupError < SrtBaseException
% SReachTools/SrtSetupError: Custom exception object for Socbox setup errors
% ============================================================================
% 
% Customized class for generating Socbox setup errors, subclass of the 
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
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 
    
    properties (Constant, Access = private)
        mnemonic = 'setupError';
    end
    
    methods
        function obj = SrtSetupError(varargin)
            obj@SrtBaseException(SrtInvalidArgsError.mnemonic, varargin{:});
       end
    end
    
    methods (Static)
        function id = getErrorId()
            id = [SrtBaseException.getErrorComponent(), ':', ...
                SrtInvalidArgsError.mnemonic];
        end
    end
end

