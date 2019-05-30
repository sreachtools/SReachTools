function bounded_set = SReachSetLagBset(sys, onestep_prob_thresh, options)
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
%        https://sreachtools.github.io/license/
%

    valid_bound_set_methods = {'load','polytope','ellipsoid'};
        % Input parsing
    valid_prob = {'term'};
    valid_method= {'chance-open','genzps-open','lag-under','lag-over'};
    validatestring(options.bound_set_method, valid_bound_set_methods);

    inpar = inputParser();
    inpar.addRequired('sys', @(x) validateattributes(x, ...
        {'LtiSystem','LtvSystem'}, {'nonempty'}));
    inpar.addRequired('onestep_prob_thresh', @(x) validateattributes(x,...
        {'numeric'}, {'scalar','>=',0,'<=',1}));
    inpar.addRequired('options', @(x) any(validatestring(x.bound_set_method,...
        valid_bound_set_methods)));
    
    try
        inpar.parse(sys, onestep_prob_thresh, options);
    catch err
        exc = SrtInvalidArgsError.withFunctionName();
        exc = exc.addCause(err);
        throwAsCaller(exc);
    end

    disturbance = sys.dist;
    
    % only need to validate attributes if not loading from file
    if ~strcmpi(options.bound_set_method, 'load')
        % validate that the disturbance is a RandomVector object
        % then ensure that the disturbance is Gaussian
        validateattributes(disturbance, {'RandomVector'}, ...
            {'nonempty'}, 'SReachSetLagBset', 'disturbance')
        switch disturbance.type
            case 'Gaussian'
                % check that the option.bound_set_method is valid
                validatestring(options.bound_set_method,...
                    {'polytope','ellipsoid'}, 'SReachSetLagBset',...
                    'options.bound_set_method for Gaussian disturbance');
            case 'UserDefined'
                % check that the option.bound_set_method is valid
                if ~strcmpi(options.bound_set_method,'polytope')
                    throwAsCaller(SrtInvalidArgsError(['Invalid options ',...
                        'provided.\nBounded disturbance set computation ', ...
                        'requires options.bound_set_method to be ', ...
                        '''polytope'' for UserDefined disturbance']));
                end
            otherwise
                throwAsCaller(SrtInvalidArgsError('Got an invalid disturbance'))
        end

        % check that onestep_prob_thresh is a value in [0,1]
        validateattributes(onestep_prob_thresh, {'numeric'},...
            {'nonempty', 'scalar','>=', 0, '<=', 1}, 'SReachSetLagBset',...
            'onestep_prob_thresh');
    end

    % check the option.bound_set_method and call appropriate sub-function
    switch(options.bound_set_method)
        case 'ellipsoid'            
            % Create the ellipsoid centered at the Gaussian mean and has a
            % shape matrix that is an appropriately scaled version of the
            % Gaussian covariance matrix
            r_squared = chi2inv(onestep_prob_thresh, disturbance.dim);
            ellipse_shape_mat = disturbance.cov() * r_squared;
            bounded_set = SReachEllipsoid(disturbance.mean(),ellipse_shape_mat);
            
        case 'polytope'
            validateattributes(options.template_polytope, {'Polyhedron'},...
                {'nonempty'}, 'SReachSetLagBset', ...
                'options.template_polytope for polytope option');
            % Use RandomVector/getProbPolyhedron to compute the bounded set as a
            % scaled version of a polytope | Returns a bounded_set that is
            % guaranteed to have probability > onestep_prob_thresh
            bounded_set = getBsetWithProb(disturbance,...
                options.template_polytope, onestep_prob_thresh,...
                options.desired_accuracy, options.verbose);

        case 'load'
            % load a predefined bounded set, primarily used for comparison with
            % previous works
            
            % variable argument should be the file location
            validateattributes(options.load_str, {'char'}, {'nonempty'},...
                'SReachSetLagBset', ...
                'options.load_str for load option');
            
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
            throwAsCaller(SrtInvalidArgsError(...
                'Invalid option.bound_set_method'));
    end
end
