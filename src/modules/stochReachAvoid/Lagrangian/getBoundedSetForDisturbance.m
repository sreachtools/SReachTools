function bounded_set = getBoundedSetForDisturbance(disturbance, ...
    horizon_length, beta, method, varargin)
% SReachTools/stochasticReachAvoid/getBoundedSetForDisturbance: Get bounded 
% disturbance set for approximation
% ============================================================================
%
% This function will get a bounded disturbance set used to compute robust
% reach avoid sets or robust effective target tubes.
%
% Usage: see examples/boundedDisturbanceSets.m
%
% ============================================================================
% 
% bounded_set = getBoundedSetForDisturbance(disturbance, ...
%     horizon_length, beta, method, varargin)
% 
% Inputs:
% -------
%   disturbance    - StochasticDisturbance object
%   horizon_length - Length of the time horizon
%   beta           - Probability threshold
%   method         - Method for computing bounded set
%   varargin       - Dependent upon method chosen, see below
%
%   Available methods:
%       'random' - Get an approximation of the ellipsoid using random
%                  direction choices; only usable for Gaussian-type
%                  disturbances; varargin must be an integer for the
%                  number of random directions to be used; e.g.
%           bounded_set = getBoundedSetForDisturbance(...
%               StochasticDisturbance('Gaussian', zeros(2,1), eye(2)), ...
%               4, ...
%               0.8, ...
%               'random', ...
%               100);
%
%       'box'    - Get an n-dimensional cuboid that satisfies the
%                  probability threshold; does not accept varargins;
%                  currenlty not implemented; e.g.
%           bounded_set = getBoundedSetForDisturbance(...
%               StochasticDisturbance('Gaussian', zeros(2,1), eye(2)), ...
%               4, ...
%               0.8, ...
%               'box');
%
%       'load'   - Load a predefined polyhedron bounding set; primarily
%                  used for comparison and repeatability testing; varargin
%                  must be a character array of the path to the file to
%                  load; mat files to be loaded must have specific design,
%                  see Notes section; when using load method all other 
%                  inputs are irrelevant; e.g.
%           bounded_set = getBoundedSetForDisturbance(...
%               [], ...
%               [], ...
%               [], ...
%               'load', ...
%               '/path/to/the/file/to/load/file.mat');
%
% Outputs:
% --------
%   bounded_set    - Polyhedron object
%
% Notes:
%   - When using the 'load' method the mat files must have only one variable
%     saved in the mat file and that variable must be a Polyhedron object.
%
% ============================================================================
%
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
%

    % only need to validate attributes if not loading from file
    if ~strcmp(method, 'load')
        % validate that the disturbance is a StochasticDisturbance object
        % then ensure that the disturbance is Gaussian
        validateattributes(disturbance, {'RandomVector'}, ...
            {'nonempty'})
        if ~strcmpi(disturbance.type, 'Gaussian')
            exc = SrtInvalidArgsError('Disturbance must be of type Gaussian');
            throwAsCaller(exc);
        end

        % check that the horizon is some nonzero integer
        validateattributes(horizon_length, {'numeric'}, ...
            {'nonzero', 'positive', 'integer'});

        % check that beta is a value in [0,1]
        validateattributes(beta, {'double'}, {'>=', 0, '<=', 1})

        % check that the method is a character or string array
        validateattributes(method, {'char', 'string'}, {'nonempty'})
    end

    % check the method and call appropriate sub-function
    switch(method)
        case 'random'
            % when choosing random direction need to specify the number of 
            % vectors to use
            validateattributes(varargin{1}, {'numeric'}, {'scalar', 'integer'});

            bounded_set = boundedEllipseByRandomVectors(...
                disturbance, ...
                horizon_length, ...
                beta, ...
                varargin{1});

        case 'box'
            if strcmpi(disturbance.type, 'gaussian')
                validateattributes(varargin{1}, {'numeric'}, ...
                    {'scalar', 'positive'});

                bounded_set = getBoundingBoxForGaussian(...
                    disturbance, ...
                    horizon_length, ...
                    beta, ...
                    varargin{1});
            end
        case 'optim-box'
            if strcmpi(disturbance.type, 'gaussian')
                validateattributes(varargin{1}, {'numeric'}, ...
                    {'nonempty'});

                bounded_set = getOptimizationBoxForGaussian(...
                    disturbance, ...
                    horizon_length, ...
                    beta, ...
                    varargin{1});
            end
        case 'optimization'
        case 'load'
            % load a predefined bounded set, primarily used for comparison with
            % previous works
            
            % variable argument should be the file location
            validateattributes(varargin{1}, {'char'}, {'nonempty'})
            
            if exist(varargin{1}, 'file') ~= 2
                exc = SrtInternalError(['Mat file to load does not ', ...
                    'exist on the path.']);
                throw(exc);
            end
            
            % load the mat file, loads as a struct
            ls = load(varargin{1});
            
            % look for the Polyhedron
            fnames = fields(ls);
            if length(fnames) > 1
                exc = SrtInternalError(['Mat file contains more than ', ...
                    'saved object. Please see Notes section of the help ', ...
                    'for details about how mat files used for loading ', ...
                    'must be structured.']);
                throw(exc);
            else
                bounded_set = ls.(fnames{1});
            end
            
            % validate that what was loaded from the mat is actually a
            % polyhedron
            validateattributes(bounded_set, {'Polyhedron'}, {'nonempty'});
            
            
        otherwise
            exc = SrtInvalidArgsError(['Invalid method provided, see ', ...
                'help for available methods'])
            throwAsCaller(exc);
            
    end

end

function bounded_set = boundedEllipseByRandomVectors(disturbance, ...
    horizon_length, beta, n_directions)
% SReachTools/getBoundedSetForDisturbance/boundedEllipseByRandomVectors: Get bounded 
% disturbance ellipse with random direction choices
% ============================================================================
%
% Get bounded disturbance set approximation as an Polyhedral overapproximation
% of an ellipse by selecting random directions on the surface of the ellipse
%
% Usage: Nested function
%
% ============================================================================
% 
% bounded_set = boundedEllipseByRandomVectors(disturbance, ...
%     horizon_length, beta, n_directions)
% 
% Inputs:
% -------
%   disturbance    - StochasticDisturbance object
%   horizon_length - Length of the time horizon
%   beta           - Probability threshold
%   n_directions   - Number or directions for the approximation
%
% Outputs:
% --------
%   bounded_set    - Polyhedron object
%
% ============================================================================
%
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
%

    % should probably provide a warning about the use of random direction for
    % 2-dimensional system
    % if disturbance.dim <= 2)
    %     warning(['For disturbances with dimension less than 2 random ', ...
    %         'there are more direct solutions for obtaining the ellipse. ', ...
    %         'Using ''lowdim'' option will provide faster and likely better ', ...
    %         'results.']);
    % end

    % compute the ellipsoid radii needed for obtaining desired probability
    r2 = chi2inv(beta^(1/horizon_length), 2);
    ellipse_rads = disturbance.parameters.covariance * r2;

    n = size(ellipse_rads,1);
    Ahalf = zeros(n_directions/2+n, n);
    bhalf = zeros(n_directions/2+n, 1);
    Ahalf(1:n,:) = eye(n)';
    bhalf(1:n) = diag(sqrt(ellipse_rads));
    for i = n+1:n_directions/2+n
        direction = 2 * randn(1,n) - 1;
        direction = direction ./ sqrt(direction * direction');
        Ahalf(i,:) = direction;
        bhalf(i) = direction * sqrt(ellipse_rads) * direction';
    end
    A = [Ahalf;-Ahalf];
    b = [bhalf;bhalf];

    % Create bounded polyhedron from A, b inequalities, i.e. Ax <= b
    bounded_set = Polyhedron(A,b);
    minVRep(bounded_set);

end

function poly = getOptimizationBoxForGaussian(disturbance, horizon_length, ...
    beta, center)
% SReachTools/getBoundedSetForDisturbance/getOptimizationBoxForGaussian: Get bounded 
% disturbance as box through solution of optimization problem
% ============================================================================
%
% Get bounded disturbance set approximation as an Polyhedral box obtained 
% through the solutions to the optimization process in 
%     [[paper]]
%
% Usage: Nested function
%
% ============================================================================
% 
% poly = getOptimizationBoxForGaussian(disturbance, horizon_length, ...
%     beta, center)
% 
% Inputs:
% -------
%   disturbance    - StochasticDisturbance object
%   horizon_length - Length of the time horizon
%   beta           - Probability threshold
%   center         - Center position of the box
%
% Outputs:
% --------
%   poly    - Polyhedron object
%
% ============================================================================
%
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
%

    center = disturbance.parameters.covariance^(-1/2)*(center - ...
        disturbance.parameters.mean);
    center = center;

    prob_threshold = beta^(1/horizon_length);
    
    perimeter_func = @(l) sum(l);
    
    l0 = ones(disturbance.dim, 1);
    
    l = fmincon(perimeter_func, l0, [], [], [], [], ...
        zeros(size(l0)), [], ...
        @(l) nonlinearOptimBoxConstraints(l, center, prob_threshold), ...
        optimoptions(@fmincon, 'Display', 'final'));

    poly = Polyhedron('lb', center-l/2, 'ub', center+l/2);
    
    poly = disturbance.parameters.covariance^(1/2) * poly + ...
        disturbance.parameters.mean;
end

function [c, ceq] = nonlinearOptimBoxConstraints(l, c, p)
% SReachTools/getOptimizationBoxForGaussian/nonlinearOptimBoxConstraints: Nonlinear
% constraints for getOptimizationBoxForGaussian
% ============================================================================
%
% Nonlinear constraints function for getOptimizationBoxForGaussian
%
% Usage: Nested function (getOptimizationBoxForGaussian)
%
% ============================================================================
% 
% [c, ceq] = nonlinearOptimBoxConstraints(l, c, p)
% 
% Inputs:
% -------
%   l - Side-length of box
%   c - Center of box
%   p - Probability value
%
% Outputs:
% --------
%   c   - Inequality constraints
%   ceq - Equality constraints
%
% ============================================================================
%
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
%

    c = p - mvncdf(c'-l'/2, c'+l'/2, zeros(size(l')), eye(length(l)));
    ceq = [];

end

function poly = getBoundingBoxForGaussian(disturbance, horizon_length, ...
    beta, err)
% SReachTools/getBoundedSetForDisturbance/getBoundingBoxForGaussian: Get bounded 
% disturbance as box through bisection
% ============================================================================
%
% Get bounded disturbance set approximation as an Polyhedral box obtained 
% through from bisecting solution
%
% Usage: Nested function
%
% ============================================================================
% 
% poly = getBoundingBoxForGaussian(disturbance, horizon_length, ...
%     beta, err)
% 
% Inputs:
% -------
%   disturbance    - StochasticDisturbance object
%   horizon_length - Length of the time horizon
%   beta           - Probability threshold
%   err            - Error threshold
%
% Outputs:
% --------
%   poly    - Polyhedron object
%
% ============================================================================
%
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
%

    MAX_ITERS = 10000;
    prob_threshold = beta^(1/horizon_length);
    a = 0;
    b = 1;
    
    mu = zeros(1, disturbance.dim);
    sigma = eye(disturbance.dim);
    
    center = zeros(1, disturbance.dim);
    dx_ones = ones(1, disturbance.dim);
    
    box = SimpleBox(center, b*dx_ones);
    p = box.computeGaussianProbability(mvncdf(box.vertices, mu, sigma));
    
    while p < prob_threshold + err
        b = 2 * b;
        box = SimpleBox(center, b*dx_ones);
        p = box.computeGaussianProbability(mvncdf(box.vertices, mu, sigma));
    end
    
    do_search = @(prob, i) ...
        (prob - prob_threshold > err || prob - prob_threshold < 0) && ...
        i < MAX_ITERS;
    iters = 0;
    while do_search(p, iters)
        b_new = b - (b - a) / 2;
        box = SimpleBox(center, b_new * dx_ones);
        p = box.computeGaussianProbability(mvncdf(box.vertices, mu, sigma));
        
        if p > prob_threshold + err
            b = b_new;
        else
            a = a + (b - a) / 2;
        end
        
        iters = iters + 1;
    end
    
    box = SimpleBox(center, b*dx_ones);
    
    poly = box.getPolyhedron();
    poly.minHRep;
    
    poly = disturbance.parameters.covariance^(1/2) * poly + ...
        disturbance.parameters.mean;
end   
