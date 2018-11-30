classdef RandomVector
% Create a random vector object
% ==========================================================================
%
% Defines a random vector. We currently support:
% 1. Gaussian    : Characterized by its mean vector and covariance matrix
% 2. UserDefined : Characterized by a random number generator
%
%
% Usage:
% ------
%
% % Define a Gaussian random variable of mean 0 and standard deviation 2:
% GaussianRV = RandomVector('Gaussian', 0, 2^2);
%
% % Define a Gaussian random vector of mean [0;2] and covariance matrix 
% % eye(2):
% GaussianRV = RandomVector('Gaussian', [0;2], eye(2));
% OR
% % GaussianRV = RandomVector.gaussian([0;2], eye(2));
%
% % Define a beta-distributed 3-dimensional random vector with parameters
% % A=B=10
% BetaRV = RandomVector('UserDefined', @(N) betarnd(10,10,[3 N]));
%   
% ==========================================================================
%
% RandomVector Properties:
% ------------------------
%   type       - Random vector type (string)
%   parameters - Random vector parameters (struct)
%                   Stores mean and covariance for Gaussian random vector
%   dim        - Random vector dimension (scalar)
%   generator  - Random variable realization generator function (function 
%                handle); should take single numeric input and return an 
%                n x p matrix---where n is the dimension of the random vector
%                and p is the number of realizations (input value)
%
% RandomVector Methods:
% ---------------------
%   RandomVector/RandomVector - Class constructor
%   mean                      - Get the mean of random vector. If the mean is
%                               not defined then an empirical mean is computed
%                               using n_particle realizations [Default: 1e4]
%   cov                       - Get the covariance of random vector. If the 
%                               mean is not defined, then an empirical mean is 
%                               computed using n_particle realizations. 
%                               [Default: 1e4]
%   pdf                       - Get the probability density function as an
%                               anonymous function, defined only for Gaussian RV
%   concat                    - Get the concatenated RV for a given time horizon
%   getRealizations           - Generate realizations of the random vector
%   getProbPolyhedron         - Get the probability of the random vector
%                               lying in a user-specified polyhedron (MPT
%                               object)
%   Static
%   ------
%   gaussian                  - Get a Gaussian random vector for a
%                               specified mean vector and covariance matrix
%   exponential               - Get an exponential random vector for a
%                               specified lambda vector
%
%   Apart from these methods, you can do disp(rv) as well as F*rv, where F
%   is an appropriately dimensioned matrix.
% 
% Notes:
% ------
% * MATLAB DEPENDENCY: Uses MATLAB's Statistics and Machine Learning Toolbox.
% 
% =========================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 
% 

    properties (SetAccess = immutable)
        % RandomVector/type
        % ==================================================================
        % 
        % String indicator of type of random vector
        %
        % Acceptable types:
        %   Gaussian
        % 
        type

        % RandomVector/parameters
        % ==================================================================
        % 
        % Struct containing random vector parameter information; will vary for
        % different random vector types. Currently only Gaussian random vector
        % are supported.
        %
        % Gaussian type:
        %   parameters.mean       - Mean vector (p x 1)
        %   parameters.covariance - Covariance matrix (p x p)
        % 
        parameters

        % RandomVector/dim
        % ==================================================================
        % 
        % Dimension of the random vector
        % 
        dim

    end

    properties (SetAccess = private)
        % RandomVector/generator
        % ==================================================================
        % 
        % Function to generate instances of the random variable
        % 
        generator
    end

    methods
        function obj = RandomVector(rv_type, varargin)
        %  Constructor for random vector class
        % ====================================================================
        %
        % Inputs:
        % -------
        %   rv_type  - Character array of random vector type; currently
        %              acceptable types:
        %                  'Gaussian'
        %                  'UserDefined'
        %   varargin - Arguments with the properties of the random vector. 
        %              Additional input argument change depending on type, see
        %              below.
        % 
        %          Gaussian:
        %               mu    - Gaussian mean
        %               sigma - Gaussian covariance matrix
        %          UserDefined:
        %               generator - Function handle to return realizations
        %                           of the random vector
        %
        % Outputs:
        % --------
        %   obj - Random vector object
        %
        % =====================================================================
        % 
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        % 
            
            % Check if the random vector type is a string 
            validatestring(rv_type, {'Gaussian', 'UserDefined'});

            switch(lower(rv_type))
                case 'gaussian'
                    obj.type = 'Gaussian'; 

                    if length(varargin) ~= 2 
                        throwAsCaller(SrtInvalidArgsError(['Gaussian random',...
                            ' vector definition needs the mean vector ', ...
                            'and covariance matrix']));
                    end

                    % Check if mean is a nonempty numeric column vector
                    validateattributes(varargin{1}, {'numeric'},...
                        {'column'},'RandomVector/RandomVector','mean');
                    % Check if covariance matrix is a nonempty square matrix
                    validateattributes(varargin{2}, {'numeric'}, ...
                       {'nonempty','square'},'RandomVector/RandomVector',...
                       'covariance');
                                   
                    % set parameters
                    obj.parameters.mean = varargin{1};
                    obj.parameters.covariance = varargin{2};
                    
                    % TODO: Throw warning as in ellipsoid
                    if ~issymmetric(obj.parameters.covariance)
                        % Compute the symmetric component of it
                        symm_cov_matrix = (obj.parameters.covariance +...
                            obj.parameters.covariance')/2;
                        % Max error element-wise
                        max_err = max(max(abs(obj.parameters.covariance -...
                            symm_cov_matrix)));
                        if max_err > eps
                            warning('SReachTools:runtime',sprintf(...
                                ['Non-symmetric covariance matrix made ',...
                                 'symmetric (max element-wise error: ',...
                                 '%1.3e)!'], max_err));
                        end
                        obj.parameters.covariance = symm_cov_matrix;
                    end
                    
                    % Check if the mean and covariance are of correct dimensions
                    if size(obj.parameters.mean, 1) ~=...
                            size(obj.parameters.covariance, 1)
                        throwAsCaller(SrtInvalidArgsError(['Mean and ',...
                            'covariance matrix have different dimensions']));
                    end
                        
                    % For some reason, -eps alone is not enough?
                    min_eig_val = min(eig(obj.parameters.covariance));
                    if  min_eig_val < -2*eps
                        throwAsCaller(SrtInvalidArgsError(['Covariance ',...
                            'matrix can not have negative eigenvalues']));
                    elseif min_eig_val <= eps
                        warning('SReachTools:runtime',['Creating a ',...
                            'Gaussian which might have a deterministic ',...
                            'component']);
                    end
                                        
                    % Update the dimension
                    obj.dim = size(obj.parameters.mean, 1);
                    
                    % Define an anonymous function using MATLAB's built-in 
                    % multivariate Gaussian random number generator
                    % We need to transpose the mean for mvnrnd
                    obj.generator = @(N) mvnrnd(obj.parameters.mean', ...
                        obj.parameters.covariance, N)';
                case 'userdefined'
                    obj.type = 'UserDefined'; 

                    % user defined variable should be called with the type and
                    % a generator function, length of varargin should be 1
                    if length(varargin) < 1
                        throwAsCaller(SrtInvalidArgsError(['Too few ', ...
                            'inputs for UserDefined random vector']));
                    elseif length(varargin) > 1
                        throwAsCaller(SrtInvalidArgsError(['Too many ', ...
                            'inputs for UserDefined random vector']));
                    end

                    % generator argument needs to be a function handle that 
                    % takes a single numeric argument, should output an n x p
                    % matrix where n in the dimension of the random vector
                    % and p is the number of points provided as an integer 
                    % argument
                    generator = varargin{1};
                    if ~isa(generator, 'function_handle')
                        throwAsCaller(SrtInvalidArgsError(['Generator ', ...
                            'for user defined variable must be of type ', ...
                            '''function_handle''']));
                    end

                    n_instances = ceil(10*rand());
                    try
                        X = generator(n_instances);
                    catch err
                        exc = SrtInvalidArgsError(['Error in calling ', ...
                            'generator function. See cause below.']);
                        exc = exc.addCause(err);
                        throw(exc);
                    end

                    if size(X, 2) ~= n_instances
                        throw(SrtInvalidArgsError(['Generator function ', ...
                            'did not return the appropriate number of ', ...
                            'realizations']));
                    end
                    % Set the dimension
                    obj.dim = size(X, 1);
                    
                    % Set the generator
                    obj.generator = generator;
                    
                    % set parameters (empty for now)
                    obj.parameters = [];
                otherwise
                    % Shouldn't even reach here!
                    throwAsCaller(SrtInvalidArgsError(sprintf(...
                        '%s-type random vectors is not supported', rv_type)));
            end
        end
        
        function disp(obj)
        % Override of MATLAB internal display
        % ====================================================================
        % 
        % Overriding of MATLAB built-in display function for the class
        %
        % ====================================================================
        % 
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        %
            
            fprintf('%d-dimensional %s random vector\n', obj.dim, obj.type);
        end
        
        function newobj = mtimes(obj, F)
        % Override of MATLAB multiplication command
        % ====================================================================
        % 
        % Inputs:
        % -------
        %   obj - RandomVector object
        %   F   - Linear transformation matrix for multiplication
        %
        % Outputs:
        % --------
        %   newobj - RandomVector object (F*obj)
        %
        % Notes:
        % ------
        % * While this function updates the generator for a UserDefined random
        %   vector, it is HIGHLY RECOMMENDED to redefine the random vector
        %   separately with an updated generator function to avoid nested
        %   generator functions
        % ====================================================================
        % 
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        %
            
            switch [class(obj), class(F)]
                case ['RandomVector','double']
                    % All ok
                case ['double', 'RandomVector']
                    % Need to switch the arguments
                    Ftemp = obj;
                    obj = F;
                    F = Ftemp;
                otherwise
                    throwAsCaller(SrtInvalidArgsError(sprintf(['Operation *',...
                       ' not defined between *%s, %s'], class(obj), class(F))));
            end
            if size(F, 2) ~= obj.dim
                throwAsCaller(SrtInvalidArgsError(['Mismatch ',...
                    'between dimensions']));
            end
            switch obj.type
                case 'Gaussian'
                    newobj=RandomVector('Gaussian', F*obj.parameters.mean,...
                        F*obj.parameters.covariance*F');
                case 'UserDefined'
                    newobj=RandomVector('UserDefined', @(N) F*obj.generator(N));
                otherwise
                    throwAsCaller(SrtInvalidArgsError(sprintf(...
                        ['Multiplication is not supported for %s-type ',...
                         'random vectors'], obj.type)));
            end
        end
        
        function newobj = concat(obj, time_horizon)
        % Create a concatenated random vector of length time_horizon
        % ====================================================================
        % 
        % Inputs:
        % -------
        %   obj           - RandomVector object (typically disturbance w_k)
        %   time_horizon  - The number of time steps of interest N
        %
        % Outputs:
        % --------
        %   newobj        - RandomVector object (for disturbance, it can be used 
        %                   as W = [w_0^\top w_1^\top ... w_{N-1}^\top])
        %
        % Notes:
        % ------
        % * We make the independent and identical assumption to obtain the
        %   concatenated random vector
        %
        % ====================================================================
        % 
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        %
            
            % Check if mean is a nonempty numeric column vector
            validateattributes(time_horizon, {'numeric'},...
                {'integer','>',0}, 'RandomVector/concat', 'time_horizon');
            
            switch obj.type
                case 'Gaussian'
                    muW = repmat(obj.parameters.mean, time_horizon, 1);
                    covW = kron(eye(time_horizon), obj.parameters.covariance);
                    newobj = RandomVector('Gaussian', muW, covW);
                case 'UserDefined'
                    % reshape and not repmat because we want independent
                    % samples
                    newGenerator = @(N) reshape(...
                        obj.generator(N*time_horizon), [], N);
                    newobj = RandomVector('UserDefined', newGenerator);
                otherwise
                    throwAsCaller(SrtInvalidArgsError(sprintf(...
                        ['Concatenation is not supported for %s-type ',...
                         'random vectors'], obj.type)));
            end
        end

        function xs = getRealizations(obj, n_realizations)
        % Generate n_realizations realizations of the random vector
        % ====================================================================
        % 
        % Inputs:
        % -------
        %   obj            - RandomVector object (typically disturbance w_k)
        %   n_realizations - Number of realizations desired
        %
        % Outputs:
        % --------
        %   xs             - Realizations matrix of dimension 
        %                    obj.dim x n_realizations. Each realization is given 
        %                    as a column vector.
        %
        % Notes:
        % ------
        % * In case of Gaussian random vectors, mvnrnd is used
        % * In case of UserDefined random vectors, the user-provided
        %   generator is used
        %
        % ====================================================================
        % 
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        %

            % Check if N is a nonempty numeric scalar
            validateattributes(n_realizations, {'numeric'},...
                {'integer','>',0}, 'RandomVector/getRealizations',...            
                'n_realizations');
            
            % Use the generator property to generate the realizations
            xs = obj.generator(n_realizations);
        end

        function pdf = pdf(obj)
        % Get the pdf of a random vector
        % ====================================================================
        % 
        % Inputs:
        % -------
        %   obj - RandomVector object
        %
        % Outputs:
        % --------
        %   pdf - Probability density function (anonymous function handle)
        %
        % Notes:
        % ------
        % * This code is not tested
        % * Only Random Vectors of type 'Gaussian' currently have a pdf.
        % * Other random vector types will throw an error
        % * The anonymous function used for the definition of obj.pdf transposes 
        %   the accepted column vector for using mvnpdf.
        % * RandomVector.pdf takes in arguments of the form N_points x
        %   random_vector_dim
        % 
        % ====================================================================
        % 
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        %

            switch obj.type
                case 'Gaussian'
                    % Define an anonymous function using MATLAB's built-in pdf 
                    % Transpose the mean for mvnpdf
                    pdf = @(x) mvnpdf(x, obj.mean()', obj.cov());
                otherwise
                    throwAsCaller(SrtInvalidArgsError(sprintf(...
                        ['No probability density function available for ',...
                         '%s-type random vectors'], obj.type)));
            end
        end

        function m = mean(obj, varargin)
        % Convenience method for accessing the (sample) mean of a random vector
        % ====================================================================
        % 
        % Inputs:
        % -------
        %   obj         - RandomVector object
        %   n_particles - [Optional] Number of particles to use in sample 
        %                 mean estimation [Default: 1e6]
        %
        %
        % Outputs:
        % --------
        %   m           - Mean of random vector
        %
        % Notes:
        % ------
        % * For Random Vectors of type 'Gaussian', we return the exact mean
        % * For Random Vectors of type 'UserDefined', we return MATLAB's
        %   estimated mean obtained from n_particles samples
        % 
        % ====================================================================
        % 
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        %

            if ~isempty(varargin) 
                n_particles = varargin{1};
                validateattributes(n_particles, {'numeric'},...
                    {'integer','>',0}, 'RandomVector/mean', 'n_particles');
                if length(varargin) > 1
                    throwAsCaller(SrtInvalidArgsError('Too many inputs'));
                end
            else
                n_particles = 1e6;
            end

            switch obj.type
                case 'Gaussian'
                    m = obj.parameters.mean;
                case 'UserDefined'
                    m = mean(obj.getRealizations(n_particles), 2);
                otherwise
                    throwAsCaller(SrtInvalidArgsError(sprintf(...
                        ['No mean available for ',...
                         '%s-type random vectors'], obj.type)));
            end
        end

        function covar = cov(obj, varargin)
        % Convenience method for accessing the covariance of a random vector
        % ====================================================================
        % 
        % Inputs:
        % -------
        %   obj         - RandomVector object
        %   n_particles - [Optional] Number of particles to use in sample 
        %                 covariance estimation [Default: 1e6]
        %
        % Outputs:
        % --------
        %   covar       - Covariance of random vector
        %
        % Notes:
        % ------
        % * For Random Vectors of type 'Gaussian', we return the exact
        %   covariance matrix
        % * For Random Vectors of type 'UserDefined', we return MATLAB's
        %   estimated covariance matrix obtained from n_particles samples
        % 
        % ====================================================================
        % 
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        %

            if ~isempty(varargin) 
                n_particles = varargin{1};
                validateattributes(n_particles, {'numeric'},...
                    {'integer','>',0}, 'RandomVector/cov', 'n_particles');
                if length(varargin) > 1
                    throwAsCaller(SrtInvalidArgsError('Too many inputs'));
                end
            else
                n_particles = 1e6;
            end

            switch obj.type
                case 'Gaussian'
                    covar = obj.parameters.covariance;
                case 'UserDefined'
                    covar = cov(obj.getRealizations(n_particles)');
                otherwise 
                    throwAsCaller(SrtInvalidArgsError(sprintf(...
                        ['No mean available for ',...
                         '%s-type random vectors'], obj.type)));
            end
        end
        
        function prob = getProbPolyhedron(obj, test_polyhedron, varargin)
        % Compute the probability of a random vector lying in the given
        % polyhedron
        % ====================================================================
        % 
        % Inputs:
        % -------
        %   obj             - RandomVector object
        %   test_polyhedron - Polyhedron object (polytope whose probability of
        %                     occurrence is of interest)
        %   n_particles     - [Optional] Number of particles to use in integration
        %                     for UserDefined random vector [Default: 1e6]
        %
        % Outputs:
        % --------
        %   covar           - Probability of the random vector lying in the
        %                     given polytope
        %
        % Notes:
        % ------
        % * For Random Vectors of type 'Gaussian', we use Genz's algorithm.
        %   We enforce an accuracy of 1e-3 via iteratedQscvmnv or a maximum of
        %   10 iterations.
        % * For Random Vectors of type 'UserDefined', we use a Monte-Carlo
        %   simulation using n_particles
        % 
        % ====================================================================
        % 
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        %

            if ~isempty(varargin) 
                n_particles = varargin{1};
                validateattributes(n_particles, {'numeric'},...
                    {'integer','>',0});
                if length(varargin) > 1
                    throwAsCaller(SrtInvalidArgsError('Too many inputs'));
                end
            else
                n_particles = 1e6;
            end

            % Check if we got a MPT's Polyhedron object of correct dimensions
            validateattributes(test_polyhedron, {'Polyhedron'}, {'nonempty'},...
                'RandomVector/getProbPolyhedron','test_polyhedron');
            if test_polyhedron.Dim ~= obj.dim
                throwAsCaller(SrtInvalidArgsError(['Mismatch in polytope ',...
                    'dimensions and random vector dimensions']));
            end

            switch obj.type
                case 'Gaussian'
                    % Use Genz's algorithm's error estimate to enforce this
                    % accuracy via iteratedQscmvnv
                    desired_accuracy = 1e-3;
                    
                    % Construct the half-space representation for qscmvnv
                    qscmvnv_lb = repmat(-Inf, [size(test_polyhedron.A, 1), 1]);
                    qscmvnv_coeff_matrix = test_polyhedron.A;
                    % We have the polytope given by Ax <= b, but qscmvnv expects 
                    % a zero mean Gaussian. Hence, define x -> x + mean
                    qscmvnv_ub = test_polyhedron.b-test_polyhedron.A*obj.mean();

                    % Call Genz's algorithm in an iterative approach to compute
                    % the probability. Uses the desired_accuracy and the
                    % error_estimate from qscmvnv to navigate the number of
                    % particles used in qscmvnv
                    try
                        prob = iteratedQscmvnv(obj.cov(), ...
                                               qscmvnv_lb, ...
                                               qscmvnv_coeff_matrix, ...
                                               qscmvnv_ub, ...
                                               desired_accuracy, ...
                                               10);
                    catch
                        %TODO-Test: Not sure when qscmvnv will bug out!
                        throw(SrtDevError(['Error in qscmvnv (Quadrature ', ...
                            'of multivariate Gaussian)']));
                    end
                case 'UserDefined'
                    mcarlo_sims = obj.getRealizations(n_particles);
                    count_contains = test_polyhedron.contains(mcarlo_sims);
                    prob = sum(count_contains)/n_particles;
                otherwise 
                    throwAsCaller(SrtInvalidArgsError(sprintf(...
                        ['Probability computation not available for ',...
                         '%s-type random vectors'], obj.type)));
            end
        end
    end

    methods (Static)
        function rv = exponential(mu)
        % Convenience method for creating exponential random vectors
        % ====================================================================
        % 
        % Static method used to conveniently create n-dimensional exponential
        % random vectors. Vectors are generated using MATLAB's exprnd function.
        % 
        % ====================================================================
        % 
        % Inputs:
        % -------
        %   mu  - Exponential mean | Must be a column vector
        %
        % Outputs:
        % --------
        %   rv - RandomVector object
        %
        % ====================================================================
        % 
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        %

            validateattributes(mu, {'numeric'},...
                    {'column','nonempty'}, 'RandomVector/exponential', 'mu');
            dim = length(mu);
            
            rv = RandomVector('UserDefined', @(N) exprnd( ...
                repmat(mu, 1, N), dim, N));
        end

        function rv = gaussian(mu, covar)
        % Convenience method for creating Gaussian random vectors
        % ====================================================================
        % 
        % Inputs:
        % -------
        %   mu     - Gaussian mean vector (column vector)
        %   covar  - Gaussian covariance matrix (square matrix)
        %
        % Outputs:
        % --------
        %   rv     - RandomVector object
        %
        % ====================================================================
        % 
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
        % 
        %

            validateattributes(mu, {'numeric'},...
                    {'column','nonempty'}, 'RandomVector/gaussian', 'mu');
            validateattributes(covar, {'numeric'},...
                    {'square','nonempty'}, 'RandomVector/gaussian', 'covar');
               
            rv = RandomVector('Gaussian', mu, covar);
        end
    end
end
