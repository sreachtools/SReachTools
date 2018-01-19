% Description : Unit test script for the LtiSystem Class
% 1/25/2018
%   - Tests all errors, asserts, and one internal error within the class (11)
%   - Tests clean LtiSystem definitions (8)
%   - Skips 1 internal errors and testing disp

% Parameter for double integrator dynamics
T = 0.5;

%% Incorrect LtiSystem definition with invalid number of arguments
correct_error_id_sent_out = 0;
try
    sys = LtiSystem('InputMatrix');
catch ME
    switch ME.identifier
        case 'SReach:invalidArgs'
            if strcmp(ME.message, ...
                      'Arguments must be given as name-value pairs')
                correct_error_id_sent_out = 1;
            else
                error('SReach:internal',...
                      'Unexpected message')
            end
        otherwise
            disp(ME)
    end
end
assert(correct_error_id_sent_out == 1,...
       'Trying an invalid number of arguments succeeded');

%% Incorrect empty LtiSystem definition is not allowed
correct_error_id_sent_out = 0;
try
    sys = LtiSystem();
catch ME
    switch ME.identifier
        case 'SReach:invalidArgs'
            if strcmp(ME.message,...
                      'State matrix can not be empty')
                correct_error_id_sent_out = 1;
            else
                error('SReach:internal',...
                      'Unexpected message')
            end
        otherwise
            disp(ME)
    end
end
assert(correct_error_id_sent_out == 1,...
       'Trying empty system definition (not allowed) succeeded');

%% Incorrect LtiSystem definition without state matrix
correct_error_id_sent_out = 0;
try
    sys = LtiSystem('InputMatrix',eye(2));
catch ME
    switch ME.identifier
        case 'SReach:invalidArgs'
            if strcmp(ME.message,...
                      'State matrix can not be empty')
                correct_error_id_sent_out = 1;
            else
                error('SReach:internal',...
                      'Unexpected message')
            end
        otherwise
            disp(ME)
    end
end
assert(correct_error_id_sent_out == 1,...
       'Trying system definition without state_matrix succeeded');

%% Incorrect LtiSystem definition with unhandled argument
correct_error_id_sent_out = 0;
try
    sys = LtiSystem('InputMatrixGoneBad',eye(2));
catch ME
    switch ME.identifier
        case 'SReach:invalidArgs'
            if strcmp(ME.message,...
                      'Unhandled argument given')
                correct_error_id_sent_out = 1;
            else
                error('SReach:internal',...
                      'Unexpected message')
            end
        otherwise
            disp(ME)
    end
end
assert(correct_error_id_sent_out == 1,...
       'Trying an unhandled argument succeeded');

%% Incorrect LtiSystem definition with invalid input matrix (wrong rows)
correct_error_id_sent_out = 0;
try
    sys = LtiSystem('StateMatrix', [1, T; 0, 1], ...
                    'InputMatrix', [T^2], ...
                    'InputSpace', Polyhedron('lb', -1, 'ub', 1));
catch ME
    switch ME.identifier
        case 'SReach:invalidArgs'
            if strcmp(ME.message,...
                      'Input matrix does not have correct row numbers')
                correct_error_id_sent_out = 1;
            else
                error('SReach:internal',...
                      'Unexpected message')
            end
        otherwise
            disp(ME)
    end
end
assert(correct_error_id_sent_out == 1,...
       'Trying an invalid input matrix (incorrect rows) succeeded');

%% Incorrect LtiSystem definition with invalid input matrix (wrong columns)
correct_error_id_sent_out = 0;
try
    sys = LtiSystem('StateMatrix', [1, T; 0, 1], ...
                    'InputMatrix', [T^2;T], ...
                    'InputSpace', Polyhedron('lb', [-1;-1], 'ub', [1;1]));
catch ME
    switch ME.identifier
        case 'SReach:invalidArgs'
            if strcmp(ME.message,...
                      'Input matrix does not have correct column numbers')
                correct_error_id_sent_out = 1;
            else
                error('SReach:internal',...
                      'Unexpected message')
            end
        otherwise
            disp(ME)
    end
end
assert(correct_error_id_sent_out == 1,...
       'Trying an invalid input matrix (incorrect columns) succeeded');

%% Incorrect LtiSystem definition with invalid input property combination
correct_error_id_sent_out = 0;
try
    sys = LtiSystem('StateMatrix', [1, T; 0, 1], ...
                    'InputMatrix', ones(2,4));
catch ME
    switch ME.identifier
        case 'SReach:invalidArgs'
            if strcmp(ME.message,...
                      'Empty input space: But non-column input matrix')
                correct_error_id_sent_out = 1;
            else
                error('SReach:internal',...
                      'Unexpected message')
            end
        otherwise
            disp(ME)
    end
end
assert(correct_error_id_sent_out == 1,...
       'Trying a empty input space with non-column input matrix succeeded');


%% Incorrect LtiSystem definition with invalid disturbance matrix (wrong rows)
correct_error_id_sent_out = 0;
try
    sys = LtiSystem('StateMatrix', [1, T; 0, 1], ...
                    'DisturbanceMatrix', [T^2], ...
                    'Disturbance', Polyhedron('lb', -1, 'ub', 1));
catch ME
    switch ME.identifier
        case 'SReach:invalidArgs'
            if strcmp(ME.message,...
                      'Disturbance matrix does not have correct row numbers')
                correct_error_id_sent_out = 1;
            else
                error('SReach:internal',...
                      'Unexpected message')
            end
        otherwise
            disp(ME)
    end
end
assert(correct_error_id_sent_out == 1,...
       'Trying an invalid disturbance matrix (incorrect rows) succeeded');

%% Incorrect LtiSystem definition with invalid disturbance matrix (wrong columns)
correct_error_id_sent_out = 0;
try
    sys = LtiSystem('StateMatrix', [1, T; 0, 1], ...
                    'DisturbanceMatrix', [T^2;T], ...
                    'Disturbance', Polyhedron('lb', [-1,-1], 'ub', [1,1]));
catch ME
    switch ME.identifier
        case 'SReach:invalidArgs'
            if strcmp(ME.message,...
                      'Disturbance matrix does not have correct column numbers')
                correct_error_id_sent_out = 1;
            else
                error('SReach:internal',...
                      'Unexpected message')
            end
        otherwise
            disp(ME)
    end
end
assert(correct_error_id_sent_out == 1,...
       'Trying an invalid disturbance matrix (incorrect columns) succeeded');

%% Incorrect LtiSystem definition with non-square state matrix
correct_error_id_sent_out = 0;
try
    sys = LtiSystem('StateMatrix', [1; 1]);
catch ME
    switch ME.identifier
        case 'SReach:invalidArgs'
            if strcmp(ME.message,...
                      'State matrix is not square')
                correct_error_id_sent_out = 1;
            else
                error('SReach:internal',...
                      'Unexpected message')
                %fprintf('Wrong error message\nExpected: ')
                %disp(ME.message)
            end
        otherwise
            disp(ME)
    end
end
assert(correct_error_id_sent_out == 1,...
       'Trying a non-square state matrix succeeded');


%% Incorrect LtiSystem definition with invalid disturbance property combination
correct_error_id_sent_out = 0;
try
    sys = LtiSystem('StateMatrix', [1, T; 0, 1], ...
                    'DisturbanceMatrix', ones(2,4));
catch ME
    switch ME.identifier
        case 'SReach:invalidArgs'
            if strcmp(ME.message,...
                      'Empty disturbance: But non-column disturbance matrix')
                correct_error_id_sent_out = 1;
            else
                error('SReach:internal',...
                      'Unexpected message')
            end
        otherwise
            disp(ME)
    end
end
assert(correct_error_id_sent_out == 1,...
       ['Trying a empty disturbance space with non-column distubance matrix',...
       ' succeeded']);

%% Incorrect LtiSystem definition with a invalid stochastic disturbance (Gauss)
correct_error_id_sent_out = 0;
mean_disturbance = zeros(5,1);
covariance_disturbance = eye(5);
GaussianDisturbance = StochasticDisturbance('Gaussian',...
                                             mean_disturbance,...
                                             covariance_disturbance);
try
    sys = LtiSystem('StateMatrix', [1, T; 0, 1], ...
                    'DisturbanceMatrix', ones(2,4), ...
                    'Disturbance', GaussianDisturbance);
catch ME
    switch ME.identifier
        case 'SReach:invalidArgs'
            if strcmp(ME.message,...
                      'Disturbance matrix does not have correct column numbers')
                correct_error_id_sent_out = 1;
            else
                error('SReach:internal',...
                      'Unexpected message')
            end
        otherwise
            disp(ME)
    end
end
assert(correct_error_id_sent_out == 1,...
       ['Trying an invalid Gaussian disturbance (dimension mismatch) succeeded']);

%% Partially Correct LtiSystem definition (Interval disturbance
%permitted without specifying input_matrix)
incorrect_error_id_sent_out = 0;
try
    sys = LtiSystem('StateMatrix', [1, T; 0, 1], ...
                    'Disturbance', Polyhedron('lb', -1, 'ub', [1]));
catch ME
    incorrect_error_id_sent_out = 1;
    error('SReach:internal',...
          'Was not expecting an error')
    throw(ME)
end
assert(incorrect_error_id_sent_out == 0,...
       'Raised error when given a valid LtiSystem definition');

%% Partially Correct LtiSystem definition (Interval input
%permitted without specifying input_matrix)
incorrect_error_id_sent_out = 0;
try
    sys = LtiSystem('StateMatrix', [1, T; 0, 1], ...
                    'InputSpace', Polyhedron('lb', -1, 'ub', [1]));
catch ME
    incorrect_error_id_sent_out = 1;
    error('SReach:internal',...
          'Was not expecting an error')
    throw(ME)
end
assert(incorrect_error_id_sent_out == 0,...
       'Raised error when given a valid LtiSystem definition');

%% Correct LtiSystem definition (only input)
incorrect_error_id_sent_out = 0;
try
    sys = LtiSystem('StateMatrix', [1, T; 0, 1], ...
                    'InputMatrix', [T^2;T], ...
                    'InputSpace', Polyhedron('lb', -1, 'ub', 1));
catch ME
    incorrect_error_id_sent_out = 1;
    error('SReach:internal',...
          'Was not expecting an error')
    throw(ME)
end
assert(incorrect_error_id_sent_out == 0,...
       'Raised error when given a valid disturbance-free LtiSystem definition');


%% Correct LtiSystem definition (only disturbance)
incorrect_error_id_sent_out = 0;
try
    sys = LtiSystem('StateMatrix', [1, T; 0, 1], ...
                    'DisturbanceMatrix', [T^2;T], ...
                    'Disturbance', Polyhedron('lb', -1, 'ub', 1));
catch ME
    incorrect_error_id_sent_out = 1;
    error('SReach:internal',...
          'Was not expecting an error')
    throw(ME)
end
assert(incorrect_error_id_sent_out == 0,...
       'Raised error when given a valid input-free LtiSystem definition');


%% Correct LtiSystem definition --- arbitrary system
incorrect_error_id_sent_out = 0;
try
    sys = LtiSystem('StateMatrix', zeros(2,2), ...
                    'InputMatrix', ones(2,4), ...
                    'InputSpace', Polyhedron('lb', -ones(4,1), 'ub', ones(4,1)),...
                    'DisturbanceMatrix', ones(2,6), ...
                    'Disturbance', Polyhedron('lb', -ones(6,1), 'ub', ones(6,1)));
catch ME
    incorrect_error_id_sent_out = 1;
    error('SReach:internal',...
          'Was not expecting an error')
    throw(ME)
end
assert(incorrect_error_id_sent_out == 0,...
       'Raised error when given a valid LtiSystem definition');


%% Correct LtiSystem definition --- double integrator
incorrect_error_id_sent_out = 0;
try
    sys = LtiSystem('StateMatrix', [1, T; 0, 1], ...
                    'InputMatrix', [T^2/2;T], ...
                    'InputSpace', Polyhedron('lb', -1, 'ub', 1),...
                    'DisturbanceMatrix', [T^2/2;T], ...
                    'Disturbance', Polyhedron('lb', -1, 'ub', 1));
catch ME
    incorrect_error_id_sent_out = 1;
    error('SReach:internal',...
          'Was not expecting an error')
    throw(ME)
end
assert(incorrect_error_id_sent_out == 0,...
       'Raised error when given a valid LtiSystem definition');

%% Correct LtiSystem definition (Gauss)
incorrect_error_id_sent_out = 0;
mean_disturbance = zeros(5,1);
covariance_disturbance = eye(5);
GaussianDisturbance = StochasticDisturbance('Gaussian',...
                                             mean_disturbance,...
                                             covariance_disturbance);
try
    sys = LtiSystem('StateMatrix', [1, T; 0, 1], ...
                    'DisturbanceMatrix', ones(2,5), ...
                    'Disturbance', GaussianDisturbance);
catch ME
    incorrect_error_id_sent_out = 1;
    error('SReach:internal',...
          'Was not expecting an error')
    throw(ME)
end
assert(incorrect_error_id_sent_out == 0,...
       'Raised error when given a valid LtiSystem definition');

%% Correct LtiSystem definition --- double integrator Gaussian perturbed
incorrect_error_id_sent_out = 0;
mean_disturbance = 0;
covariance_disturbance = 1;
GaussianDisturbance = StochasticDisturbance('Gaussian',...
                                             mean_disturbance,...
                                             covariance_disturbance);
try
    sys = LtiSystem('StateMatrix', [1, T; 0, 1], ...
                    'InputMatrix', [T^2/2;T], ...
                    'InputSpace', Polyhedron('lb', -1, 'ub', 1),...
                    'DisturbanceMatrix', [T^2/2;T], ...
                    'Disturbance', GaussianDisturbance);
catch ME
    incorrect_error_id_sent_out = 1;
    error('SReach:internal',...
          'Was not expecting an error')
    throw(ME)
end
assert(incorrect_error_id_sent_out == 0,...
       'Raised error when given a valid LtiSystem definition');
