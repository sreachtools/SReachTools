% Description : Unit test script for the StochasticDisturbance Class
% 1/26/2018
%   - Tests all errors and asserts within the class (2 of them)
%   - Tests 2 clean StochasticDisturbance definitions
%   - Skips 1 internal errors and testing disp(object)

%% Incorrect StochasticDisturbance definition incomplete Gaussian definition
correct_error_id_sent_out = 0;
try
    GaussianDisturbance = StochasticDisturbance('Gaussian', ...
                                                 zeros(5,1));
catch ME
    switch ME.identifier
        case 'SReachTools:invalidArgs'
            if strcmp(ME.message, ...
                      ['Gaussian disturbance needs the mean vector and ', ...
                      'covariance matrix'])
                correct_error_id_sent_out = 1;
            else
                disp(ME.message)
                error('SReachTools:internal', ...
                      'Unexpected message')
            end
        otherwise
            disp(ME)
    end
end
assert(correct_error_id_sent_out == 1, ...
       'Trying an invalid number of arguments succeeded');

%% Incorrect StochasticDisturbance definition Gauss parameter dimension mismatch
correct_error_id_sent_out = 0;
try
    GaussianDisturbance = StochasticDisturbance('Gaussian', ...
                                                 zeros(5,1), ...
                                                 eye(4));
catch ME
    switch ME.identifier
        case 'SReachTools:invalidArgs'
            if strcmp(ME.message, ...
                      'Mean and covariance matrix have different dimensions')
                correct_error_id_sent_out = 1;
            else
                disp(ME.message)
                error('SReachTools:internal', ...
                      'Unexpected message')
            end
        otherwise
            disp(ME)
    end
end
assert(correct_error_id_sent_out == 1, ...
       'Trying an invalid number of arguments succeeded');

%% Correct StochasticDisturbance definition (1-dim Gaussian, match pdfs)
incorrect_error_id_sent_out = 0;
try
    GaussianDisturbance = StochasticDisturbance('Gaussian', ...
                                                 1, ...
                                                 4);
    assert( all(abs(normpdf([-10:0.1:10]',1,2) - ...
                            GaussianDisturbance.pdf([-10:0.1:10]'))<1e-8), ...
            'Pdfs mismatch was not expected');
    assert(GaussianDisturbance.dim == 1, ...
            'Dimension mismatch between the object and the covariance matrix');

catch ME
    disp(ME.message)
    incorrect_error_id_sent_out = 1;
    error('SReachTools:internal', ...
          'Was not expecting an error')
    throw(ME)
end
assert(incorrect_error_id_sent_out == 0, ...
       'Raised error when given a valid StochasticDisturbance definition');

%% Correct StochasticDisturbance definition (n-dim Gaussian, match pdfs)
incorrect_error_id_sent_out = 0;
mean_disturbance = zeros(5,1);
covariance_disturbance = eye(5);
sample_points = randn(100,5);
try
    GaussianDisturbance = StochasticDisturbance('Gaussian', ...
                                                 mean_disturbance, ...
                                                 covariance_disturbance);
    assert( all(abs(mvnpdf(sample_points, ...
                           mean_disturbance', ...
                           covariance_disturbance) - ...
                      GaussianDisturbance.pdf(sample_points))<1e-8...
               ), ...
            'Pdfs mismatch was not expected');
    assert(GaussianDisturbance.dim == size(covariance_disturbance, 2), ...
            'Dimension mismatch between the object and the covariance matrix');
catch ME
    incorrect_error_id_sent_out = 1;
    error('SReachTools:internal', ...
          'Was not expecting an error')
    throw(ME)
end
assert(incorrect_error_id_sent_out == 0, ...
       'Raised error when given a valid StochasticDisturbance definition');
