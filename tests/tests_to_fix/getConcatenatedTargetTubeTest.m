% Description : Unit test script for the function
%               @LtiSystem/getConcatTargetTube
% 1/27/2018
%   - Tests four asserts, positive scalar time_horizon and non-empty
%     polyhedron for safe set and target set, dimension match (2+2+2+1)

%% Incorrect non-scalar time horizon given
correct_error_id_sent_out = 0;
safe_set = Polyhedron('lb', -1, 'ub', 1);
target_set = Polyhedron('lb', -1, 'ub', 1);
try
    [concat_target_tube_A, concat_target_tube_b] = ...
                            getConcatTargetTube(safe_set, ...
                                                      target_set, ...
                                                      [1,1])
catch ME
    switch ME.identifier
        case 'SReachTools:invalidArgs'
            if strcmp(ME.message, ...
                      'Expected a scalar positive time_horizon')
                correct_error_id_sent_out = 1;
            else
                error('SReachTools:internal', ...
                      'Unexpected message')
            end
        otherwise
            disp(ME)
    end
end
assert(correct_error_id_sent_out == 1, ...
       'Trying a non-scalar time horizon succeeded');

%% Incorrect zero time horizon given
correct_error_id_sent_out = 0;
safe_set = Polyhedron('lb', -1, 'ub', 1);
target_set = Polyhedron('lb', -1, 'ub', 1);
try
    [concat_target_tube_A, concat_target_tube_b] = ...
                            getConcatTargetTube(safe_set, ...
                                                      target_set, ...
                                                      0)
catch ME
    switch ME.identifier
        case 'SReachTools:invalidArgs'
            if strcmp(ME.message, ...
                      'Expected a scalar positive time_horizon')
                correct_error_id_sent_out = 1;
            else
                error('SReachTools:internal', ...
                      'Unexpected message')
            end
        otherwise
            disp(ME)
    end
end
assert(correct_error_id_sent_out == 1, ...
       'Trying zero time horizon succeeded');

%% Incorrect empty safe set
correct_error_id_sent_out = 0;
safe_set = Polyhedron();
target_set = Polyhedron('lb', -1, 'ub', 1);
try
    [concat_target_tube_A, concat_target_tube_b] = ...
                            getConcatTargetTube(safe_set, ...
                                                      target_set, ...
                                                      10)
catch ME
    switch ME.identifier
        case 'SReachTools:invalidArgs'
            if strcmp(ME.message, ...
                      'Safe set must be a non-empty polyhedron')
                correct_error_id_sent_out = 1;
            else
                error('SReachTools:internal', ...
                      'Unexpected message')
            end
        otherwise
            disp(ME)
    end
end
assert(correct_error_id_sent_out == 1, ...
       'Trying empty safe set succeeded');

%% Incorrect empty target set
correct_error_id_sent_out = 0;
safe_set = Polyhedron('lb', -1, 'ub', 1);
target_set = Polyhedron();
try
    [concat_target_tube_A, concat_target_tube_b] = ...
                            getConcatTargetTube(safe_set, ...
                                                      target_set, ...
                                                      10)
catch ME
    switch ME.identifier
        case 'SReachTools:invalidArgs'
            if strcmp(ME.message, ...
                      'Target set must be a non-empty polyhedron')
                correct_error_id_sent_out = 1;
            else
                error('SReachTools:internal', ...
                      'Unexpected message')
            end
        otherwise
            disp(ME)
    end
end
assert(correct_error_id_sent_out == 1, ...
       'Trying empty target set succeeded');

%% Incorrect diff dimensions target and safe set
correct_error_id_sent_out = 0;
safe_set = Polyhedron('lb', -1, 'ub', 1);
target_set = Polyhedron('lb', [-1,-1], 'ub', [1,1]);
try
    [concat_target_tube_A, concat_target_tube_b] = ...
                            getConcatTargetTube(safe_set, ...
                                                      target_set, ...
                                                      10)
catch ME
    switch ME.identifier
        case 'SReachTools:invalidArgs'
            if strcmp(ME.message, ...
                      'Safe and target sets must be of the same dimension');
                correct_error_id_sent_out = 1;
            else
                error('SReachTools:internal', ...
                      'Unexpected message')
            end
        otherwise
            disp(ME)
    end
end
assert(correct_error_id_sent_out == 1, ...
       'Trying different dimensional safe and target set succeeded');

%% Incorrect safe set type
correct_error_id_sent_out = 0;
safe_set = StochasticDisturbance('Gaussian',0,1);
target_set = Polyhedron('lb', [-1,-1], 'ub', [1,1]);
try
    [concat_target_tube_A, concat_target_tube_b] = ...
                            getConcatTargetTube(safe_set, ...
                                                      target_set, ...
                                                      10)
catch ME
    switch ME.identifier
        case 'SReachTools:invalidArgs'
            if strcmp(ME.message, ...
                      'Safe set must be a non-empty polyhedron')
                correct_error_id_sent_out = 1;
            else
                error('SReachTools:internal', ...
                      'Unexpected message')
            end
        otherwise
            disp(ME)
    end
end
assert(correct_error_id_sent_out == 1, ...
       'Trying a wrong object for safe set succeeded');

%% Incorrect target set type
correct_error_id_sent_out = 0;
target_set = StochasticDisturbance('Gaussian',0,1);
safe_set = Polyhedron('lb', [-1,-1], 'ub', [1,1]);
try
    [concat_target_tube_A, concat_target_tube_b] = ...
                            getConcatTargetTube(safe_set, ...
                                                      target_set, ...
                                                      10)
catch ME
    switch ME.identifier
        case 'SReachTools:invalidArgs'
            if strcmp(ME.message, ...
                      'Target set must be a non-empty polyhedron')
                correct_error_id_sent_out = 1;
            else
                error('SReachTools:internal', ...
                      'Unexpected message')
            end
        otherwise
            disp(ME)
    end
end
assert(correct_error_id_sent_out == 1, ...
       'Trying a wrong object for target set succeeded');
