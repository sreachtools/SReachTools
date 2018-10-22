function setSrtWarning(warn_str,state)
    valid_warn_str = {'SReachTools:runtime','SReachTools:desiredAccuracy', ...
        'SReachTools:setup'};
    if strcmpi(warn_str,'all') && length(state) == length(valid_warn_str)
        for warn_indx = 1:length(valid_warn_str)
            warn_str = valid_warn_str{warn_indx};
            warning(state{warn_indx},warn_str);
        end
    elseif strcmpi(warn_str,'all') && any(contains({'on','off'},state)) == 1
        for warn_indx = 1:length(valid_warn_str)
            warn_str = valid_warn_str{warn_indx};
            warning(state,warn_str);
        end
    elseif any(contains(valid_warn_str,warn_str))...
            && any(contains({'on','off'},state)) == 1
        warning(state,warn_str);
    else
        throw(SrtInvalidArgsError('Invalid warning string(s)/states provided'));
    end
end