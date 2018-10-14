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
% bounded_set = getBoundedSetForDisturbance(disturbance, time_horizon,...
%   onestep_prob_thresh, option)
%                                        e.g.
%           bounded_set = getBoundedSetForDisturbance(...
%               RandomVector('Gaussian', zeros(2,1), eye(2)), ...
%               4, ...
%               0.8, ...
%               'random', ...
%               100);
%           bounded_set = getBoundedSetForDisturbance(...
%               RandomVector('Gaussian', zeros(2,1), eye(2)), ...
%               4, ...
%               0.8, ...
%               'box');
%           bounded_set = getBoundedSetForDisturbance(...
%               [], ...
%               [], ...
%               [], ...
%               'load', ...
%               '/path/to/the/file/to/load/file.mat');
%
% 
% Inputs:
% -------
%   disturbance - RandomVector object
%   time_horizon- Length of the time horizon
%   onestep_prob_thresh - Probability threshold
%   option      - TODO
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

    valid_bound_set_methods = {'load','random','box','optim-box'};
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
        validateattributes(options.bound_set_method, {'char', 'string'},...
            {'nonempty'})
    end

    % check the option.bound_set_method and call appropriate sub-function
    switch(options.bound_set_method)
        case 'random'
            % when choosing random direction need to specify the number of 
            % vectors to use
            validateattributes(options.num_dirs, {'numeric'},...
                {'scalar', 'integer'});

            bounded_set = boundedEllipseByRandomVectors(disturbance,...
                onestep_prob_thresh, options.num_dirs);

        case 'box'
            % Has to be Gaussian
            validateattributes(options.err_thresh, {'numeric'}, ...
                {'scalar', 'positive'});

            bounded_set = getBoundingBoxForGaussian(disturbance, ...
                onestep_prob_thresh, options.err_thresh);
        case 'optim-box'
            % Has to be Gaussian
            validateattributes(options.box_center, {'numeric'}, {'nonempty',...
                'vector','size', [disturbance.dim 1]});
            
            bounded_set = getOptimizationBoxForGaussian(disturbance, ...
                onestep_prob_thresh, options.box_center);
        case 'load'
            % load a predefined bounded set, primarily used for comparison with
            % previous works
            
            % variable argument should be the file location
            validateattributes(options.load_str, {'char'}, {'nonempty'})
            
            if exist(options.load_str, 'file') ~= 2
                throwAsCaller(SrtInvalidArgsError(['Mat file to load does ',...
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
%                 throw(SrtInvalidArgsError(['Mat file did not contain a ',...
%                     'polyhedron of the desired type']));
%             end            
        otherwise
            exc = SrtInvalidArgsError('Invalid option.bound_set_method');
            throwAsCaller(exc);            
    end

end

function bounded_set = boxBisection(disturbance, onestep_onestep_prob_thresh, err)

    MAX_ITERS = 10000;
    % xi2_onestep_prob_thresh = onestep_prob_thresh^(1/time_horizon);
    a = 0;
    b = 1;
    
    mu = zeros(1, disturbance.dim);
    sigma = eye(disturbance.dim);
    
    center = zeros(1, disturbance.dim);
    dx_ones = ones(1, disturbance.dim);

    bx = Polyhedron('A', [eye(size(c, 1)), -eye(size(c, 1))], ...
                    'b', b * ones(2 * size(c, 1), 1));
    bx = disturbance.parameters.covariance^(1/2) * bx + ...
        disturbance.parameters.mean;
    p = computeProb(disturbance, bx);
    
    while p < onestep_onestep_prob_thresh + err
        b = 2 * b;
        bx = Polyhedron('A', [eye(size(c, 1)), -eye(size(c, 1))], ...
                        'b', b * ones(2 * size(c, 1), 1));
        bx = disturbance.parameters.covariance^(1/2) * bx + ...
            disturbance.parameters.mean;
        p = computeProb(disturbance, bx);
    end
    
    do_search = @(prob, i) ...
        (prob - onestep_onestep_prob_thresh > err || prob - onestep_onestep_prob_thresh < 0) && ...
        i < MAX_ITERS;
    iters = 0;
    while (p - onestep_onestep_prob_thresh > err || ...
           p - onestep_onestep_prob_thresh < 0 ) && ...
           iters < MAX_ITERS

        % bisect 
        b_new = b - (b - a) / 2;
        bx = Polyhedron('A', [eye(size(c, 1)), -eye(size(c, 1))], ...
                        'b', b * ones(2 * size(c, 1), 1));
        bx = disturbance.parameters.covariance^(1/2) * bx + ...
            disturbance.parameters.mean;
        p = computeProb(disturbance, bx);
        
        if p > onestep_onestep_prob_thresh + err
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
%  Get bounded 
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
    %     warning(['For disturbances with dimension less than 2 random ', ...
    %         'there are more direct solutions for obtaining the ellipse. ', ...
    %         'Using ''lowdim'' option will provide faster and likely ',...
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

function poly = getOptimizationBoxForGaussian(disturbance, ...
    onestep_prob_thresh, center)
%  Get bounded 
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
% poly = getOptimizationBoxForGaussian(disturbance, time_horizon, ...
%     onestep_prob_thresh, center)
% 
% Inputs:
% -------
%   disturbance    - RandomVector object
%   time_horizon - Length of the time horizon
%   onestep_prob_thresh           - Probability threshold
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
    
    xi2_onestep_prob_thresh = onestep_prob_thresh;
    
    perimeter_func = @(l) sum(l);
    
    l0 = ones(disturbance.dim, 1);
    
    l = fmincon(perimeter_func, l0, [], [], [], [], ...
        zeros(size(l0)), [], ...
        @(l) nonlinearOptimBoxConstraints(l, center, xi2_onestep_prob_thresh), ...
        optimoptions(@fmincon, 'Display', 'final'));

    poly = Polyhedron('lb', center-l/2, 'ub', center+l/2);
    
    poly = disturbance.parameters.covariance^(1/2) * poly + ...
        disturbance.parameters.mean;
end

function [c, ceq] = nonlinearOptimBoxConstraints(l, c, p)
%  Nonlinear
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

function poly = getBoundingBoxForGaussian(disturbance, ...
    onestep_prob_thresh, err)
%  Get bounded 
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
% poly = getBoundingBoxForGaussian(disturbance, time_horizon, ...
%     onestep_prob_thresh, err)
% 
% Inputs:
% -------
%   disturbance    - RandomVector object
%   time_horizon - Length of the time horizon
%   onestep_prob_thresh           - Probability threshold
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
    xi2_onestep_prob_thresh = onestep_prob_thresh;
    a = 0;
    b = 1;
    
    mu = zeros(1, disturbance.dim);
    sigma = eye(disturbance.dim);
    
    center = zeros(1, disturbance.dim);
    dx_ones = ones(1, disturbance.dim);
    
    bx = boxFromCenterAndHalfLength(center, b);
    bx = disturbance.parameters.covariance^(1/2) * bx + ...
        disturbance.parameters.mean;
    p = computeProb(disturbance, bx);
    
    while p < xi2_onestep_prob_thresh + err
        b = 2 * b;
        box = SimpleBox(center, b*dx_ones);
        p = box.computeGaussianProbability(mvncdf(box.vertices, mu, sigma));
    end
    
    do_search = @(prob, i) ...
        (prob - xi2_onestep_prob_thresh > err || prob - xi2_onestep_prob_thresh < 0) && ...
        i < MAX_ITERS;
    iters = 0;
    while do_search(p, iters)
        b_new = b - (b - a) / 2;
        box = SimpleBox(center, b_new * dx_ones);
        p = box.computeGaussianProbability(mvncdf(box.vertices, mu, sigma));
        
        if p > xi2_onestep_prob_thresh + err
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


function prob = computeProb(dist, bset)

    % Construct the half-space representation for qscmvnv
    cov_mat = (dist.parameters.covariance +...
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

function b = boxFromCenterAndHalfLength(c, hl)

    v = zeros(2^size(c, 1), size(c, 1));

    if size(c, 1) == 2
        P = [ hl(1),  hl(2); ...
             -hl(1),  hl(2); ...
              hl(1), -hl(2); ...
             -hl(1), -hl(2)];

        v = c + P;
    elseif size(c, 1) == 3
        P = [ hl(1),  hl(2),  hl(3); ...
              hl(1), -hl(2),  hl(3); ...
              hl(1),  hl(2), -hl(3); ...
              hl(1), -hl(2), -hl(3); ...
             -hl(1),  hl(2),  hl(3); ...
             -hl(1), -hl(2),  hl(3); ...
             -hl(1),  hl(2), -hl(3); ...
             -hl(1), -hl(2), -hl(3)];

        v = c + P;
    elseif size(c, 1) == 4
        P = [ hl(1),  hl(2),  hl(3),  hl(4); ...
              hl(1),  hl(2),  hl(3), -hl(4); ...
              hl(1),  hl(2), -hl(3),  hl(4); ...
              hl(1),  hl(2), -hl(3), -hl(4); ...
              hl(1), -hl(2),  hl(3),  hl(4); ...
              hl(1), -hl(2),  hl(3), -hl(4); ...
              hl(1), -hl(2), -hl(3),  hl(4); ...
              hl(1), -hl(2), -hl(3), -hl(4); ...
             -hl(1),  hl(2),  hl(3),  hl(4); ...
             -hl(1),  hl(2),  hl(3), -hl(4); ...
             -hl(1),  hl(2), -hl(3),  hl(4); ...
             -hl(1),  hl(2), -hl(3), -hl(4); ...
             -hl(1), -hl(2),  hl(3),  hl(4); ...
             -hl(1), -hl(2),  hl(3), -hl(4); ...
             -hl(1), -hl(2), -hl(3),  hl(4); ...
             -hl(1), -hl(2), -hl(3), -hl(4)];

        v = c + P;
    else
        sv = -ones(1, size(c, 1));
        lv = 1;
        while true
            done = false;
            v(lv, :) = c + hl .* p;

            for rv = size(c, 1):-1:1
                if sv(rv) < -1;
                    sv(rv) == 1;

                    if rv == 1
                        done = true;
                    else
                        sv(rv - 1) = sv(rv - 1) - 1;
                    end
                end
            end

            if done
                break;
            end
        end
    end

    b = Polyhedron(v);
    b = minHRep();
end