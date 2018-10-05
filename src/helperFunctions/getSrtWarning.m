function state = getSrtWarning()
    struct_temp = warning('query', 'SReachTools:runtime');
    state{1} = struct_temp.state;
    struct_temp = warning('query', 'SReachTools:desiredAccuracy');
    state{2} = struct_temp.state;
end