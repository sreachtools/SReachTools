% Description : Unit test script for the function
%               @LtiSystem/getConcatInputSpace
% 1/27/2018
%   - Tests two asserts, positive scalar time_horizon and non-empty
%     polyhedron for input_space (3)
%   - Tests one clean cases and compares it with the expected results (1)

umax = 1;
%% Incorrect empty polyhedral input space
sys = LtiSystem('StateMatrix', eye(2));
correct_error_id_sent_out = 0;
try
    [concatenated_input_space_A, concatenated_input_space_b] = ...
                                         getConcatInputSpace(sys,...
                                                                   2);
catch ME
    switch ME.identifier
        case 'SReachTools:invalidArgs'
            if strcmp(ME.message, ...
                      'Expected a non-empty polyhedral input space')
                correct_error_id_sent_out = 1;
            else
                error('SReachTools:internal',...
                      'Unexpected message')
            end
        otherwise
            disp(ME)
    end
end
assert(correct_error_id_sent_out == 1,...
       'Trying a non-scalar time horizon succeeded');

%% Incorrect non-scalar time horizon given
sys = LtiSystem(...
    'StateMatrix', eye(2), ...
    'InputMatrix', ones(2,1), ...
    'InputSpace', Polyhedron('lb', -umax, 'ub', umax));
correct_error_id_sent_out = 0;
try
    [concatenated_input_space_A, concatenated_input_space_b] = ...
                                         getConcatInputSpace(sys,...
                                                                   [2,2]);
catch ME
    switch ME.identifier
        case 'SReachTools:invalidArgs'
            if strcmp(ME.message, ...
                      'Expected a scalar positive time_horizon')
                correct_error_id_sent_out = 1;
            else
                error('SReachTools:internal',...
                      'Unexpected message')
            end
        otherwise
            disp(ME)
    end
end
assert(correct_error_id_sent_out == 1,...
       'Trying a non-scalar time horizon succeeded');

%% Incorrect zero time horizon given
sys = LtiSystem(...
    'StateMatrix', eye(2), ...
    'InputMatrix', ones(2,1), ...
    'InputSpace', Polyhedron('lb', -umax, 'ub', umax));
correct_error_id_sent_out = 0;
try
    [concatenated_input_space_A, concatenated_input_space_b] = ...
                                         getConcatInputSpace(sys,...
                                                                   0);
catch ME
    switch ME.identifier
        case 'SReachTools:invalidArgs'
            if strcmp(ME.message, ...
                      'Expected a scalar positive time_horizon')
                correct_error_id_sent_out = 1;
            else
                error('SReachTools:internal',...
                      'Unexpected message')
            end
        otherwise
            disp(ME)
    end
end
assert(correct_error_id_sent_out == 1,...
       'Trying zero time horizon succeeded');

%% Correct case (U^N compared)
sys = LtiSystem(...
    'StateMatrix', eye(2), ...
    'InputMatrix', ones(2,1), ...
    'InputSpace', Polyhedron('lb', -umax, 'ub', umax));
time_horizon = 10;
incorrect_error_id_sent_out = 0;
try
    [concatenated_input_space_A, concatenated_input_space_b] = ...
                                         getConcatInputSpace(sys,...
                                                                   time_horizon);
    assert( Polyhedron('H',[concatenated_input_space_A,...
                            concatenated_input_space_b]) == ...
            Polyhedron('lb', -umax * ones(time_horizon,1),...
                       'ub',  umax * ones(time_horizon,1)),...
            'SReachTools:internal',...
            'Mismatch in U^N');
catch ME
    disp(ME.message)
    incorrect_error_id_sent_out = 1;
    error('SReachTools:internal',...
          'Was not expecting an error')
    throw(ME)
end
assert(incorrect_error_id_sent_out == 0,...
       ['Raised error when given a valid disturbance-free double integrator',...
       'system']);

