function state = getSrtWarning(varargin)
    valid_warn_str = {'SReachTools:runtime','SReachTools:desiredAccuracy', ...
        'SReachTools:setup'};
    if isempty(varargin)
        state = cell(length(valid_warn_str),1);
        for warn_indx = 1:length(valid_warn_str)
            warn_str = valid_warn_str{warn_indx};
            struct_temp = warning('query', warn_str);
            state{warn_indx} = struct_temp.state;        
        end
    elseif length(varargin) == 1 && any(contains(valid_warn_str,varargin{1}))
        struct_temp = warning('query', varargin{1});
        state = struct_temp.state;                
    else
        throw(SrtInvalidArgsError('Invalid warning string provided'));
    end
end