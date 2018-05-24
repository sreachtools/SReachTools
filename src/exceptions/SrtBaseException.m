classdef SrtBaseException < MException
% SReachTools/SrtBaseException: Custom exception object for Socbox internal
% errors
% ============================================================================
% 
% Customized class for generating Socbox internal errors, subclass of the 
% standard MATLAB MException class
%
% Usage:
% ------
% exc = SrtBaseException('error message')
%
% ============================================================================
%
% See also MException
%
% ============================================================================
%
%   This function is part of the Stochastic Optimal Control Toolbox.
%   License for the use of this function is given in
%        https://github.com/abyvinod/SReachTools/blob/master/LICENSE
% 
    
    properties (Access = private)
        mnemonic = '';
    end
    
    methods
        function obj = SrtBaseException(id, varargin)
            obj@MException(sprintf('SReachTools:%s', id), ''); 
            if length(varargin) >= 1
                obj.message = sprintf(varargin{:});
            end

            obj.mnemonic = id;
        end
    end

    methods (Static)
        function comp = getErrorComponent()
            comp = 'SReachTools';
        end
    end
end

