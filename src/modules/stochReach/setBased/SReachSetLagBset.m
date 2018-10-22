function bounded_set = SReachSetLagBset(disturbance, onestep_prob_thresh, ...
    options)
% Get bounded disturbance set for approximation
% ============================================================================
%
% This function will get a bounded disturbance set used to compute robust
% reach avoid sets or robust effective target tubes.
%
% Usage: see examples/boundedDisturbanceSets.m
%
% ============================================================================
% 
% bounded_set = SReachSetLagBset(disturbance, time_horizon, ...
%   onestep_prob_thresh, option)
% 
% Inputs:
% -------
%   disturbance         - RandomVector object
%   time_horizon        - Length of the time horizon
%   onestep_prob_thresh - Probability threshold
%   option              - Struct specifying methods to obtain the bounded set
%                         see SReachSetOptions
%
% Outputs:
% --------
%   bounded_set    - Polyhedron object
%
%   
%   See also SReachSetOptions
% 
% Notes:
% * When using the 'load' method the mat files must have only one variable
%   saved in the mat file and that variable must be a Polyhedron object.
%
% ============================================================================
%
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
%

    valid_bound_set_methods = {'load','random','box'};
    validatestring(options.bound_set_method, valid_bound_set_methods);

    % only need to validate attributes if not loading from file
    if ~strcmpi(options.bound_set_method, 'load')
        % validate that the disturbance is a RandomVector object
        % then ensure that the disturbance is Gaussian
        validateattributes(disturbance, {'RandomVector'}, ...
            {'nonempty'})
        if ~strcmpi(disturbance.type, 'Gaussian')
            exc = SrtInvalidArgsError('Disturbance must be of type Gaussian');
            throwAsCaller(exc);
        end

        % check that the horizon is some nonzero integer
        % validateattributes(time_horizon, {'numeric'}, ...
        %     {'nonzero', 'positive', 'integer'});

        % check that onestep_prob_thresh is a value in [0,1]
        validateattributes(onestep_prob_thresh, {'double'}, {'>=', 0, '<=', 1})

        % check that the option.bound_set_method is a character or string array
        validateattributes(options.bound_set_method, {'char', 'string'}, ...
            {'nonempty'})
    end

    % check the option.bound_set_method and call appropriate sub-function
    switch(options.bound_set_method)
        case 'random'
            % when choosing random direction need to specify the number of 
            % vectors to use
            validateattributes(options.num_dirs, {'numeric'}, ...
                {'scalar', 'integer'});

            bounded_set = boundedEllipseByRandomVectors(disturbance, ...
                onestep_prob_thresh, options.num_dirs);

        case 'box'
            % Has to be Gaussian
            validateattributes(options.err_thresh, {'numeric'}, ...
                {'scalar', 'positive'});

            bounded_set = boxBisection(disturbance, ...
                onestep_prob_thresh, options.err_thresh);
        case 'optim-box'
            % Has to be Gaussian
            validateattributes(options.box_center, {'numeric'}, {'nonempty', ...
                'vector','size', [disturbance.dim 1]});
            
            bounded_set = getOptimizationBoxForGaussian(disturbance, ...
                onestep_prob_thresh, options.box_center);
        case 'load'
            % load a predefined bounded set, primarily used for comparison with
            % previous works
            
            % variable argument should be the file location
            validateattributes(options.load_str, {'char'}, {'nonempty'})
            
            if exist(options.load_str, 'file') ~= 2
                throwAsCaller(SrtInvalidArgsError(['Mat file to load does ', ...
                    'not exist on the path.']));
            end
            
            % load the mat file, loads as a struct
            ls = load(options.load_str);
            
            % look for the Polyhedron
            fnames = fields(ls);
            if length(fnames) > 1
                exc = SrtInvalidArgsError(['Mat file contains more than ', ...
                    'saved object. Please see Notes section of the help ', ...
                    'for details about how mat files used for loading ', ...
                    'must be structured.']);
                throwAsCaller(exc);
            else
                bounded_set = ls.(fnames{1});
            end
            
            % validate that what was loaded from the mat is actually a
            % polyhedron
            validateattributes(bounded_set, {'Polyhedron'}, {'nonempty'});
%           TODO-Test: Need to check for disturbance dimension            
%             if bounded_set.Dim ~= disturbance.dim
%                 throw(SrtInvalidArgsError(['Mat file did not contain a ', ...
%                     'polyhedron of the desired type']));
%             end            
        otherwise
            exc = SrtInvalidArgsError('Invalid option.bound_set_method');
            throwAsCaller(exc);            
    end

end

function bounded_set = boxBisection(disturbance, onestep_prob_threshold, err)
% Get bounded set that is a n-d box via bisection method
% =============================================================================
% 
% Get box bounded disturbance set for Lagrangian approximation using bisection
% method. Box will be centered at the mean of the disturbance and have
% its shape altered by the covariance of the disturbance. For i.i.d. 
% disturbances this will create an n-d rectangular object. The polytope will 
% have 2 * n facets and 2^n vertices
% 
% Usage: Nested function of SReachSetLagBSet
% 
% =============================================================================
% 
% bounded_set = boxBisection(disturbance, onestep_prob_threshold, err)
% 
% Inputs:
% -------
%   disturbance            - RandomVector object
%   onestep_prob_threshold - One-step probability threshold. For 
%                            underapproximation:
%                               level_set_threshold^(1 / time_horizon)
%                            For overapproximation
%                               (1 - level_set_threshold)^(1 / time_horizon)
%   err - Error tolerance
% 
% Outputs:
% --------
%   bounded_set - Polyhedron object
% 
% ============================================================================
%
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
%


    MAX_ITERS = 10000;
    % xi2_onestep_prob_thresh = onestep_prob_thresh^(1/time_horizon);
    a = 0;
    b = 1;
    
    mu = zeros(1, disturbance.dim);
    sigma = eye(disturbance.dim);
    
    center = zeros(1, disturbance.dim);
    dx_ones = ones(1, disturbance.dim);

    bx = Polyhedron('A', [eye(disturbance.dim); -eye(disturbance.dim)], ...
                    'b', b * ones(2 * disturbance.dim, 1));
    bx = disturbance.parameters.covariance^(1/2) * bx + ...
        disturbance.parameters.mean;
    p = computeProb(disturbance, bx);
    
    while p < onestep_prob_threshold + err
        b = 2 * b;
        bx = Polyhedron('A', [eye(disturbance.dim); -eye(disturbance.dim)], ...
                        'b', b * ones(2 * disturbance.dim, 1));
        bx = disturbance.parameters.covariance^(1/2) * bx + ...
            disturbance.parameters.mean;
        p = computeProb(disturbance, bx);
    end
    
    do_search = @(prob, i) ...
        (prob - onestep_prob_threshold > err || ...
         prob - onestep_prob_threshold < 0) && ...
        i < MAX_ITERS;
    iters = 0;
    while (p - onestep_prob_threshold > err || ...
           p - onestep_prob_threshold < 0 ) && ...
           iters < MAX_ITERS

        % bisect 
        b_new = b - (b - a) / 2;
        bx = Polyhedron('A', [eye(disturbance.dim); -eye(disturbance.dim)], ...
                        'b', b_new * ones(2 * disturbance.dim, 1));
        bx = disturbance.parameters.covariance^(1/2) * bx + ...
            disturbance.parameters.mean;
        p = computeProb(disturbance, bx);
        
        if p > onestep_prob_threshold + err
            b = b_new;
        else
            a = a + (b - a) / 2;
        end
        
        iters = iters + 1;
    end
    
    bounded_set = bx;

end

function bounded_set = boundedEllipseByRandomVectors(disturbance, ...
    onestep_prob_thresh, n_directions)
% Get bounded disturbance ellipse with random direction choices
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
%     time_horizon, onestep_prob_thresh, n_directions)
% 
% Inputs:
% -------
%   disturbance    - RandomVector object
%   time_horizon - Length of the time horizon
%   onestep_prob_thresh           - Probability threshold
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
    %     warn me(['For disturbances with dimension less than 2 random ', ...
    %         'there are more direct solutions for obtaining the ellipse. ', ...
    %         'Using ''lowdim'' option will provide faster and likely ', ...
    %         'better results.']);
    % end

    % compute the ellipsoid radii needed for obtaining desired probability
    r2 = chi2inv(onestep_prob_thresh, 2);
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

function prob = computeProb(dist, bset)
% Compute probability that distrubance lies in polyhedral object
% ============================================================================
%
% Compute probability that distrubance lies in polyhedral object
%
% Usage: Nested function
%
% ============================================================================
% 
% prob = computeProb(dist, bset)
% 
% Inputs:
% -------
%   dist - RandomVector object
%   bset - Bounded set (Polyhedron object)
%
% Outputs:
% --------
%   prob - Probablity that dist lies in bset
%
% ============================================================================
%
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
%

    % Construct the half-space representation for qscmvnv
    cov_mat = (dist.parameters.covariance + ...
        dist.parameters.covariance')/2; 
    qscmvnv_lb = repmat(-Inf, [size(bset.A, 1), 1]);
    qscmvnv_coeff_matrix = bset.A;
    qscmvnv_ub = bset.b - bset.A * dist.parameters.mean;
    prob = iteratedQscmvnv(cov_mat, ...
                           qscmvnv_lb, ...
                           qscmvnv_coeff_matrix, ...
                           qscmvnv_ub, ...
                           1e-3, ...
                           10);            
end