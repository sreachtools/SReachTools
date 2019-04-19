classdef SrtInvalidArgsError < SrtBaseException
% Custom exception object for SReachTools invalid arguments errors
% ============================================================================
% 
% Customized class for generating SReachTools invalid arguments errors, subclass 
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
%        https://sreachtools.github.io/license/
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
        % Throw invalid args and provide function name that received the 
        % invalid call
        % =====================================================================
        % 
        % Method to preformat error strings to throw to the user the specific
        % function that was used when passing the invalid args error. Will 
        % display message:
        % 
        %       Invalid arguments provided to << function name >>
        % 
        % =====================================================================
        % 
            stk = dbstack(1);
            try
                stk = stk(1);
            catch err
                warning('SReachTools:runtime',['Dbstack empty, perhaps ', ...
                    'you are calling SrtInvalidArgsError from command ', ...
                    'window?']);
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

