classdef SrtBaseException < MException
% Custom exception object for SReachTools internal errors
% ============================================================================
% 
% Customized class for generating SReachTools internal errors, subclass of the 
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
%        https://sreachtools.github.io/license/
% 
    
    properties (Access = private)
        % SrtBaseException/mnemonic
        % ====================================================================
        % 
        % MATLAB error mnemonic, used for building error ids
        % 
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

