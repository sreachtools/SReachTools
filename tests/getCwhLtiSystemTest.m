% Description : Unit test script for the function
%               systemDefinitions/getCwhLtiSystem()
% 1/27/2018
%   - Tests the input handling for the function (5 tests)
%   - Tests a valid input configuration (1 test)

GaussianDist = StochasticDisturbance('Gaussian',...
                                     zeros(4,1),...
                                     diag([1e-4, 1e-4, 5e-8, 5e-8]));
inputSpace = Polyhedron('lb', -0.01*ones(2,1),...
                        'ub',  0.01*ones(2,1));                                 
%% Incorrect invalid input space (non-polytope)
correct_error_id_sent_out = 0;
try
   sys = getCwhLtiSystem(4,... 
                         'wer',...
                         GaussianDist);
catch ME
    switch ME.identifier
        case 'SReachTools:invalidArgs'
            if strcmp(ME.message, ...
                      'Must provide polyhedral input space')
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
       'Trying an invalid input space (non-polytope) succeeded');

%% Incorrect invalid input space (incorrect dimensional polytope)
correct_error_id_sent_out = 0;
try
   sys = getCwhLtiSystem(4, ...
                         Polyhedron('lb', -0.01*ones(3,1),...
                                    'ub',  0.01*ones(3,1)),...
                         GaussianDist);
catch ME
    switch ME.identifier
        case 'SReachTools:invalidArgs'
            if strcmp(ME.message, ...
                      'Input matrix does not have correct column numbers')
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
       ['Trying an invalid input space (incorrect dimensional  polytope)',...
       ' succeeded']);

%% Incorrect dimensioned Gaussian Disturbance provided
% Intricacies of stochatic disturbance tested within itself
correct_error_id_sent_out = 0;
try
   sys = getCwhLtiSystem(4,...
                         inputSpace,...
                         StochasticDisturbance(...
                                    'Gaussian',...
                                    zeros(2,1),...
                                    diag([1e-4, 1e-4])));
catch ME
    switch ME.identifier
        case 'SReachTools:invalidArgs'
            if strcmp(ME.message, ...
                      'Disturbance matrix does not have correct column numbers')
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
       'Trying an invalid mean vector (incorrect column vector) succeeded');


%% Correct function call
incorrect_error_id_sent_out = 0;
try
   sys = getCwhLtiSystem(4,...
                         inputSpace,...
                         GaussianDist);
catch ME
    disp(ME.message)
    incorrect_error_id_sent_out = 1;
    error('SReachTools:internal',...
          'Was not expecting an error')
    throw(ME)
end
assert(incorrect_error_id_sent_out == 0,...
       'Raised error when given a correct function call');

%% TODO: Untested 6 dim, user params