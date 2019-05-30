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
% % OR
% GaussianRV = RandomVector.gaussian([0;2], eye(2));
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
%   Apart from these methods, for a RandomVector object rv, you can do:
%   disp(rv)                  - Display information about rv
%   F*rv, rv*F                - Both work as F x rv (NO TRANSPOSE enforced)
%                               where F is an appropriately dimensioned matrix
%   rv + v                    - Add v, a deterministic vector of appropriate
%                               dimension, to rv
%   [rv1;rv2;rv3]             - Concatenate multiple random vectors
% 
% Notes:
% ------
% * MATLAB DEPENDENCY: Uses MATLAB's Statistics and Machine Learning Toolbox.
%
% 
% =========================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://sreachtools.github.io/license/
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
        %      https://sreachtools.github.io/license/
        % 
        % 

            % Tolerances
            % Maximum elementwise error tolerance when making a matrix symmetric
            max_elem_err_tol = 1e-10;           
            % Saturation for minimum eigenvalues
            min_eig_val_saturation = 1e-10;
            
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
                    
                    % Ensure that the covariance matrix is symmetric
                    if ~issymmetric(obj.parameters.covariance)
                        % Compute the symmetric component of it
                        symm_cov_matrix = (obj.parameters.covariance +...
                            obj.parameters.covariance')/2;
                        % Max error element-wise
                        max_err = max(max(abs(obj.parameters.covariance -...
                            symm_cov_matrix)));
                        if max_err > max_elem_err_tol
                            throwAsCaller(SrtInvalidArgsError(sprintf(...
                                ['Non-symmetric covariance matrix provided ',...
                                 '(max element-wise error: %1.3e)'], max_err)));
                        elseif max_err > eps
                            warning('SReachTools:runtime',sprintf(...
                                ['Non-symmetric covariance matrix made ',...
                                 'symmetric (max element-wise error: ',...
                                 '%1.3e)!'], max_err));
                        end
                        obj.parameters.covariance = symm_cov_matrix;
                    end
                    
                    % Ensure that the covariance matrix has real
                    % nonnegative eigenvalues
                    % For some reason, -eps alone is not enough?
                    % TODO: We can do better using eigs | However, there is a
                    %       compatibility issue for MATLAB earlier than 2017a
                    % min_eig_val = eigs(obj.parameters.covariance,1,'smallestreal');
                    min_eig_val = min(eig(obj.parameters.covariance));
                    if min_eig_val <= -2*eps
                        throwAsCaller(SrtInvalidArgsError(['Covariance ',...
                            'matrix can not have negative eigenvalues']));
                    elseif min_eig_val < 0
                        warning('SReachTools:runtime', sprintf(['Sanitized ', ...
                            'covariance matrix since negative eigenvalues ',...
                            '> -2*eps and <0 found!\nNew covariance matrix ',...
                            'has all the eigenvalues below %1.0e set to 0.'],...
                            min_eig_val_saturation));
                        [V, E] = eig(obj.parameters.covariance);
                        eig_vector = diag(E);
                        eig_vector_sanitized = zeros(size(eig_vector));
                        nnz_eig_vector=(abs(eig_vector)>min_eig_val_saturation);
                        eig_vector_sanitized(nnz_eig_vector) = ...
                            eig_vector(nnz_eig_vector);
                        E_sanitized = diag(eig_vector_sanitized);
                        cov_temp = V * E_sanitized * V';                     
                        % Sometimes symmetricity is lost. So just to be
                        % sure.
                        obj.parameters.covariance = (cov_temp + cov_temp')/2;
                    end
                    
                    

                    % Check if the mean and covariance are of correct dimensions
                    if size(obj.parameters.mean, 1) ~=...
                            size(obj.parameters.covariance, 1)
                        throwAsCaller(SrtInvalidArgsError(['Mean and ',...
                            'covariance matrix have different dimensions']));
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
        %      https://sreachtools.github.io/license/
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
        %      https://sreachtools.github.io/license/
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
                       ' not defined between %s and %s'], class(obj),...
                       class(F))));
            end
            if isequal(size(F),[1 1])
                % F is a scalar
                F = F * eye(obj.dim);
            end
            if size(F, 2) ~= obj.dim
                throwAsCaller(SrtInvalidArgsError(['Mismatch ',...
                    'between dimensions']));
            end
            switch obj.type
                case 'Gaussian'
                    newobj=RandomVector('Gaussian',F*obj.mean(),F*obj.cov()*F');
                case 'UserDefined'
                    newobj=RandomVector('UserDefined', @(N) F*obj.generator(N));
                otherwise
                    throwAsCaller(SrtInvalidArgsError(sprintf(...
                        ['Multiplication is not supported for %s-type ',...
                         'random vectors'], obj.type)));
            end
        end
        
        function newobj = plus(obj, v)
        % Override of MATLAB plus command
        % ====================================================================
        % 
        % Inputs:
        % -------
        %   obj - RandomVector object
        %   v   - Deterministic vector to be added to the random vector OR
        %         a RandomVector object
        %
        % Outputs:
        % --------
        %   newobj - RandomVector object (obj + v)
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
        %      https://sreachtools.github.io/license/
        % 
        %
            
            summands_type = [];
            switch [class(obj), class(v)]
                case ['RandomVector','double']
                    % Check dimensions
                    if ~isequal(size(v), [obj.dim 1])
                        throwAsCaller(SrtInvalidArgsError(['Mismatch in ',...
                            'dimensions of the random vector and v']));
                    end
                    % Set the flag for the type of summation
                    summands_type = 'determ_vec_plus_rv';
                case ['double', 'RandomVector']
                    % Need to switch the arguments
                    vtemp = obj;
                    obj = v;
                    v = vtemp;
                    % Check dimensions
                    if ~isequal(size(v), [obj.dim 1])
                        throwAsCaller(SrtInvalidArgsError(['Mismatch in ',...
                            'dimensions of the random vector and v']));
                    end
                    % Set the flag for the type of summation
                    summands_type = 'determ_vec_plus_rv';
                case ['RandomVector', 'RandomVector']
                    % Check dimensions
                    if v.dim ~= obj.dim
                        throwAsCaller(SrtInvalidArgsError(['Mismatch in ',...
                            'dimensions of the random vector and v']));
                    end
                    % Set the flag for the type of summation
                    summands_type = 'rv_plus_rv';
                otherwise
                    throwAsCaller(SrtInvalidArgsError(sprintf(['Operation +',...
                       ' not defined between %s and %s'], class(obj),...
                       class(v))));
            end
            switch summands_type
                case 'determ_vec_plus_rv'
                    if isequal(size(v),[1 1])
                        % F is a scalar
                        v = v * ones(obj.dim,1);
                    end
                    switch obj.type
                        case 'Gaussian'
                            newobj=RandomVector.gaussian(v + obj.mean(),...
                                obj.cov());
                        case 'UserDefined'
                            newobj=RandomVector('UserDefined',...
                                @(N) repmat(v, 1, N) + obj.generator(N));
                        otherwise
                            throwAsCaller(SrtInvalidArgsError(sprintf(...
                                ['Plus is not supported for %s-type ',...
                                 'random vectors'], obj.type)));
                    end
                case 'rv_plus_rv'
                    switch [obj.type, v.type]
                        case ['Gaussian', 'Gaussian']
                            newobj=RandomVector.gaussian(...
                                v.mean() + obj.mean(), v.cov() + obj.cov());
                        case {['Gaussian', 'UserDefined'], ...
                                ['UserDefined','Gaussian'], ...
                                ['UserDefined','UserDefined']}
                            newobj=RandomVector('UserDefined',...
                                @(N) obj.generator(N) + v.generator(N));
                        otherwise
                            throwAsCaller(SrtInvalidArgsError(sprintf(...
                                ['Plus is not supported for %s-type and %s-',...
                                 'type random vectors'], obj.type, v.type)));
                    end
                otherwise
                    % Will never come here
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
        %      https://sreachtools.github.io/license/
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
        %      https://sreachtools.github.io/license/
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
        %      https://sreachtools.github.io/license/
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
        %      https://sreachtools.github.io/license/
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
        % Compute the probability of a random vector lying in a polyhedron 
        % ====================================================================
        % The probability computation is done via Monte-Carlo (quasi Monte-Carlo 
        % in Gaussian case), which is a guaranteed underapproximation (with a 
        % probabilistic guarantee). The guarantee in the general case comes
        % via Hoeffding's inequality, and, in the Gaussian case, via the 
        % confidence interval.
        % 
        % General case: A distribution-independent lower bound on the number
        % of particles needed is obtained via Hoeffding's inequality.
        %
        % For \beta = exp(-2 * n_particles * \delta^2), Hoeffding's inequality 
        % states that
        %
        %               Prob{X - E[X] \geq \delta} \leq \beta, 
        %
        % where, X is the empirical average of a collection of n_particles 
        % random variables bounded between [0,1], E[X] is the true mean, \delta 
        % is the desired_accuracy, \beta is the failure risk (the probability of 
        % the statement "X - E[X] \geq \delta" fails. 
        %
        % In this function, we set the bounded random variables to be Bernoulli 
        % random variables that is 1 when a realization of the random vector 
        % is within the given test_polyhedron, and 0 otherwise.
        % Consequently, E[X] simplifies to the probability of the
        % realization to lie in the polyhedron. 
        %
        % Given a failure risk \beta and desired accuracy \delta, we backcompute 
        % the number of particles as the following
        % 
        %               n_particles = -ln(\beta) / (2 * delta^2) 
        % 
        % This enforces the condition that empirical probability estimate X does 
        % not overapproximate the true probability E[X] by more than delta.
        %
        % Gaussian case: We use Genz's algorithm to estimate the
        % probability. We increase the number of particles in the powers of
        % 10, till error estimate is within the desired tolerance.
        %
        % In both of the cases, to remain conservative in our estimate, we 
        % subtract the desired_accuracy from the emperical probability estimate.
        % In other words, we use the lower bound of the confidence interval.
        % Thus, the probability estimate is guaranteed to be an
        % underapproximation.
        % 
        % Inputs:
        % -------
        %   obj             - RandomVector object
        %   test_polyhedron - Polyhedron object (polytope whose probability of
        %                     occurrence is of interest)
        %   desired_accuracy- [Optional] Maximum absolute deviation from the
        %                     true probability estimate [Default: 1e-2]
        %
        % Outputs:
        % --------
        %   covar           - Probability of the random vector lying in the
        %                     given polytope
        % Notes:
        % ------
        % * Due to the inverse-square dependence on the desired_accuracy, we
        %   impose a hard lower bound of 1e-2 on the desired_accuracy. This
        %   leads to the requirement of 2e5 particles.
        % * We set the default failure risk as 2e-15.
        % * Ill-formed (mean/cov has Inf/NaN) Gaussian random vectors return
        %   zero probability.
        % * We compute the floor of the probability value based on the
        %   desired_accuracy to obtain a preferable lower bound on the
        %   probability.
        %
        % ====================================================================
        % 
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://sreachtools.github.io/license/
        % 
        %

            failure_risk = 2e-15;
            gauss_scale_factor = 10;
            if ~isempty(varargin) 
                desired_accuracy = varargin{1};
                validateattributes(desired_accuracy, {'numeric'},...
                    {'scalar','>=', 1e-2, '<=', 1});
                if length(varargin) > 1
                    throwAsCaller(SrtInvalidArgsError('Too many inputs'));
                end
            else
                desired_accuracy = 1e-2;
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
                    if any(isinf(abs(obj.mean()))) || ...
                            any(isnan(abs(obj.mean()))) || ...
                            any(any(isinf(abs(obj.cov())))) || ...
                            any(any(isnan(abs(obj.cov()))))
                        % Return zero probability if ill-formed random vector
                        temp_probability = 0;
                    else
                        % Construct the half-space representation for qscmvnv
                        qscmvnv_lb = repmat(-Inf, ...
                            [size(test_polyhedron.A, 1),1]);
                        qscmvnv_coeff_matrix = test_polyhedron.A;
                        % We have the polytope given by Ax <= b, but qscmvnv 
                        % expects a zero mean Gaussian. Hence, define x->x+mean
                        qscmvnv_ub = test_polyhedron.b - ...
                            test_polyhedron.A*obj.mean();
                        % Set up for an iterative call of Genz's particles;
                        % Increase particles till the error estimate is within 
                        % threshold
                        n_particles = 10;
                        est_3sig_CI = 1;
                        while est_3sig_CI > desired_accuracy
                            try                            
                                % Call Genz's algorithm; 
                                [temp_probability, est_3sig_CI] = qscmvnv( ...
                                    n_particles,obj.cov(), qscmvnv_lb, ...
                                    qscmvnv_coeff_matrix, qscmvnv_ub);
                            catch ME
                                disp(ME);
                                throw(SrtDevError(['Error in qscmvnv (', ...
                                    'Quadrature of multivariate Gaussian)']));
                            end
                            n_particles = n_particles * gauss_scale_factor;
                        end
                    end
                case 'UserDefined'
                    % By Hoeffding's inequality
                    n_particles = ceil(-log(failure_risk)/ ...
                        (2*desired_accuracy^2));

                    if n_particles > 2e5
                        warning('SReachTools:runtime', sprintf(['Number of ',...
                                'particles required: %1.2e. Consider ', ...
                                'increasing (relaxing) desired_accuracy!'], ...
                                n_particles));
                    end

                    mcarlo_sims = obj.getRealizations(n_particles);
                    count_contains = test_polyhedron.contains(mcarlo_sims);
                    temp_probability = sum(count_contains)/n_particles;
                otherwise 
                    throwAsCaller(SrtInvalidArgsError(sprintf(...
                        ['Probability computation not available for ',...
                         '%s-type random vectors'], obj.type)));
            end
            % Subtract the desired_accuracy to remain an underapproximation
            temp_probability = temp_probability - desired_accuracy;
            
            if temp_probability > 0                 
                % Rounding DOWN the integral to the desired accuracy
                prob = floor(temp_probability/desired_accuracy) * ...
                    desired_accuracy;        
            else
                prob = 0;
            end
        end
        
        function newobj = vertcat(varargin)
        % Vertical concatentation routine
        % ====================================================================
        % 
        % Inputs:
        % -------
        %   obj             - A collection of RandomVector objects
        %
        % Outputs:
        % --------
        %   newobj          - A RandomVector object that is the vertical
        %                     concatenation of random vectors
        %
        % Notes:
        % ------
        % * This function requires the objects that are being concatenated
        %   be of same RandomVector.type.
        % * If the concatenation of a deterministic vector dv is desired with
        %   a RandomVector object rv, please use the following affine
        %   transformation-based command:
        %
        %   new_rv = [dv; zeros(rv.dim,1)] + [zeros(length(dv), rv.dim);
        %                                     eye(rv.dim) ] * rv;
        % 
        % ====================================================================
        % 
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://sreachtools.github.io/license/
        % 
    
            if nargin == 0
                throwAsCaller(SrtInvalidArgsError('Empty input provided!'));
            else
                rv_first = varargin{1};
                rv_type = rv_first.type;
                % Confirm all of them are:
                % 1. RandomVector objects
                % 2. same type
                for indx = 1:nargin
                    if ~isa(varargin{indx},'RandomVector') || ...
                            ~strcmpi(varargin{indx}.type, rv_type)
                        throwAsCaller(SrtInvalidArgsError(['Concatenation', ...
                            ' only permitted for RandomVector objects of ', ...
                            'same type.']));
                    end
                end
                switch rv_type
                    case 'Gaussian'
                        muW = [];
                        covW = [];
                        for indx = 1:nargin
                            muW = [muW;
                                   varargin{indx}.mean()];
                            new_cov = varargin{indx}.cov();
                            covW = [covW, zeros(size(covW,1), size(new_cov,2));
                                    zeros(size(new_cov,1), size(covW,2)),...
                                        new_cov];                                
                        end
                        newobj = RandomVector('Gaussian', muW, covW);
                    case 'UserDefined'
                        % Following 
                        % https://stackoverflow.com/questions/11232323/
                        % combining-anonymous-functions-in-matlab
                        % Example given there is:
                        % ca = {@(X) X, @(X) X+1, @(X) X^2};
                        % h=@(x) cellfun(@(y) y(x), ca);
                        newGeneratorCell = cell(nargin,1);
                        for indx = 1:nargin
                            newGeneratorCell{indx} = varargin{indx}.generator;
                        end
                        newGenerator = @(N) cell2mat(cellfun(@(fun) fun(N), ...
                            newGeneratorCell, 'UniformOutput', false));
                        newobj = RandomVector('UserDefined', newGenerator);
                    otherwise
                        throwAsCaller(SrtInvalidArgsError(sprintf(...
                            ['Concatenation is not supported for %s-type ',...
                             'random vectors'], obj.type)));
                end
                
            end
        end
        
        function newobj = concat(obj, repeat_times)
        % Create a concatenated random vector of length time_horizon
        % ====================================================================
        % 
        % Inputs:
        % -------
        %   obj           - RandomVector object (typically disturbance w_k)
        %   repeat_times  - The number of time steps of interest N
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
        %      https://sreachtools.github.io/license/
        % 
        %
            
            % Check if repeat_times is a nonempty positive integer
            validateattributes(repeat_times, {'numeric'},...
                {'integer','>',0}, 'RandomVector/concat', 'time_horizon');
            newobj = obj;
            for indx = 2:repeat_times
                newobj = [newobj;obj];
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
        %      https://sreachtools.github.io/license/
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
        %      https://sreachtools.github.io/license/
        % 
        %

            validateattributes(mu, {'numeric'},...
                    {'column','nonempty'}, 'RandomVector/gaussian', 'mu');
            validateattributes(covar, {'numeric'},...
                    {'square','nonempty'}, 'RandomVector/gaussian', 'covar');
               
            rv = RandomVector('Gaussian', mu, covar);
        end
        
        function rv = uniform(lb, ub)
            
            validateattributes(lb, {'numeric'}, {'column', 'nonempty'}, ...
                'RandomVector/uniform', 'lb');
            validateattributes(lb, {'numeric'}, {'column', 'nonempty'}, ...
                'RandomVector/uniform', 'ub');
            
            if any(size(lb) ~= size(ub))
                throw(SrtInvalidArgsError(['Lower and upper bounds must ', ...
                    'be the same size.']));
            end
                
            rv = RandomVector('UserDefined', @(N) ub * rand(1, N) + lb);
        end
        
        function horzcat()
        % Horizontal concatenation prevention routine!
        % ====================================================================
        % 
        % Notes:
        % ------
        % * This function just throws an error since horizontal
        %   concatenation produces random matrix!
        % 
        % ====================================================================
        % 
        % This function is part of the Stochastic Reachability Toolbox.
        % License for the use of this function is given in
        %      https://sreachtools.github.io/license/
        % 
        %
            throwAsCaller(SrtInvalidArgsError(['Horizontal concatentation ',...
                'of random vectors is not permitted!']));
        end        
    end
end
