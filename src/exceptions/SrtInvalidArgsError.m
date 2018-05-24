classdef SrtInvalidArgsError < SrtBaseException
% SReachTools/SrtInvalidArgsError: Custom exception object for Socbox 
% invalid arguments errors
% ============================================================================
% 
% Customized class for generating Socbox invalid arguments errors, subclass 
% of the standard SrtBaseException class
%
% Usage:
% ------
% exc = SrtInvalidArgsError('error message')
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
    
    properties (Constant, Access = private)
        mnemonic = 'invalidArgs';
    end
    
    methods
        function obj = SrtInvalidArgsError(varargin)
            obj@SrtBaseException(SrtInvalidArgsError.mnemonic, varargin{:});
       end
    end
    
    methods (Static)
        function id = getErrorId()
            id = [SrtBaseException.getErrorComponent(), ':', ...
                SrtInvalidArgsError.mnemonic];
        end

        function obj = withFunctionName(varargin)
            stk = dbstack(1);
            try
                stk = stk(1);
            catch err
                warning(['Dbstack empty, perhaps you are calling from ', ...
                    'the command window?']);
                stk = struct('name', '<< empty stack >>');
            end

            err_msg = sprintf('Invalid arguments provided to %s.', ...
                stk.name);

            if nargin > 0            
                err_msg = [err_msg, '\n\n', sprintf(varargin{:})];
            end

            obj = SrtInvalidArgsError(err_msg);
        end
    end
end

