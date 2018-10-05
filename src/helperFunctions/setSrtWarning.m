function setSrtWarning(state)
    if isa(state, 'cell') && length(state) > 1
        warning('SReachTools:runtime',state{1});
        warning('SReachTools:desiredAccuracy',state{2});
    else
        warning('SReachTools:runtime',state);
        warning('SReachTools:desiredAccuracy',state);
    end
end