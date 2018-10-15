function options = SReachSetOptions(prob_str, method_str, varargin)
% Create user-specifiable options for use with SReachSet()
% =============================================================================
%
% SReachSetOptions creates a MATLAB struct that contains user-specifiable
% options that may be used with SReachSet
%
% =============================================================================
%
%   options = SReachSetOptions(prob_str, method_str, varargin)
% 
% Inputs:
% -------
%   prob_str    - String specifying the problem of interest. For each case, we
%                 compute the optimal value function that maps initial states
%                 to different maximal reach probabilities
%                     1. 'first' : Stay within the safety_tube and reach the
%                                  target set early if possible
%                     2. 'term' : Stay within the safety_tube
%   method_str  - Solution technique to be used (user-specifiable
%                 options associated with each technique is enumerated)
%                     'chance-open':
%                          Convex chance-constrained approach for an open-loop
%                          controller synthesis
%                          1. pwa_accuracy: Accuracy of the piecewise affine
%                               approximation of norminvcdf used
%                     'genzps-open':
%                          Genz's algorithm + Patternsearch
%                          1. desired_accuracy: Accuracy of Gaussian integral =>
%                               Accuracy of the result
%                          2. PSoptions: MATLAB struct from psoptimset()
%                     'lag-over'/'lag-under':
%                          Lagrangian-based over- and underapproximation
%                          bound_set_method:
%                          1a. 'random'- Get an approximation of the ellipsoid
%                                        using random direction choices; only
%                                        usable for Gaussian-type disturbances;
%                                        varargin must be an integer for the
%                                        number of random directions to be used;
%                                        b. num_dirs --- Number of directions to
%                                           sample the ellipsoid for a polytopic
%                                           representation
%                          2a. 'box'    - Get an n-dimensional cuboid centered at
%                                        the disturbance mean that satisfies the
%                                        probability threshold
%                                        b. err_thresh --- Tolerance for
%                                           the bisection algorithm that
%                                           identifies the length of the box
%                          3a. 'optim-box' 
%                                      - Get an n-dimensional cuboid centered at
%                                        a user-specified box center that
%                                        satisfies the probability threshold
%                                        b. box_center --- Center for the box
%                          4a. 'load'   - Load a predefined polyhedron bounding
%                                        set; primarily used for comparison and
%                                        repeatability testing
%                                        b. load_str --- Path to the file to
%                                           load. All other inputs are
%                                           IRRELEVANT for this option.
%                                           Mat files to be loaded must
%                                           have specific design TODO
%
% Outputs:
% --------
%   options     - Collection of user-specified options for 'chance-affine'
%                 (Matlab struct created using SReachSetOptions)
%
% See also SReachSet.
%
% Notes:
% * SReachSet() will call this function internally using the default
%     values if SReachSetOptions()-based options is not explicitly provided
%     to SReachSet().
% ============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 
%

    valid_prob = {'term'};
    valid_method= {'chance-open','genzps-open','lag-under','lag-over'};
    valid_bound_method = {'load','random','box'};
    % Input parsing
    inpar = inputParser();
    inpar.addRequired('prob_str', @(x) any(validatestring(x,valid_prob)));
    inpar.addRequired('method_str', @(x) any(validatestring(x,valid_method)));

    if contains(method_str,'lag-')
        inpar.addParameter('bound_set_method',[], @(x) any(validatestring(x,...
            valid_bound_method)));
        %% Optional arguments that are made REQUIRED based on bound_set_method
        % Get number of directions for random option
        inpar.addParameter('num_dirs', 10, @(x) validateattributes(x,...
            {'numeric'}, {'scalar', 'integer'}));
        % Get bisection error threshold for box option
        inpar.addParameter('err_thresh', 1e-3, @(x) validateattributes(x,...
            {'numeric'}, {'scalar', 'positive'}));
        % Get box center for optim-box option
        inpar.addParameter('box_center', zeros(2,1),...
            @(x) validateattributes(x,{'numeric'}, {'nonempty', 'vector'}));
        % Get load string for load option
        inpar.addParameter('load_str', ' ', @(x) validateattributes(x,...
            {'char'}, {'nonempty'}));
    else
%         % TODO: Required
%         init_safe_set_affine_const, ...
%         set_of_dir_vecs,...
                                           
        switch lower(method_str)
            case 'genzps-open'            
                inpar.addParameter('desired_accuracy',1e-3, @(x)...
                    validateattributes(x, {'numeric'}, {'scalar','>',0}));
                inpar.addParameter('PSoptions',psoptimset('display','off'));

                % Ensure that patternsearch is installed
                v = ver;            
                has_fmincon = any(strcmp(cellstr(char(v.Name)),...
                    'Global Optimization Toolbox'));
                if ~has_fmincon
                    exc = SrtSetupError(['SReachSet with ''genzps-open'' ',...
                        'option needs MATLAB''s Global Optimization Toolbox.']);
                    throw(exc);
                end
            case 'chance-open'
                % Accuracy of piecewise-affine approximation of norminvcdf
                inpar.addParameter('pwa_accuracy',1e-3, @(x)...
                    validateattributes(x, {'numeric'}, {'scalar','>',0}));
                inpar.addParameter('init_safe_set_affine',Polyhedron(), @(x)...
                    validateattributes(x, {'Polyhedron'}, {'nonempty'}));
                inpar.addParameter('set_of_dir_vecs',[], @(x)...
                    validateattributes(x, {'numeric'}, {'nonempty'}));
                inpar.addParameter('verbose',0, @(x)...
                    validateattributes(x, {'numeric'},...
                    {'scalar','>=',0,'<=',1}));
        end
    end
    
    %% Construct the options MATLAB struct
    inpar.parse(prob_str, method_str, varargin{:});
    options = inpar.Results;
    
    %% Ensure that user provided an option for the must use case
    if contains(method_str,'lag-') 
        if any(contains(inpar.UsingDefaults,'bound_set_method'))
            throw(SrtInvalidArgsError(['bound_set_method is a required ',...
                'input for SReachSet when using ''lag-over''/''lag-under''.']));
        end
        switch(lower(options.bound_set_method))
            % case 'optim-box'
            %     v = ver;            
            %     has_fmincon = any(strcmp(cellstr(char(v.Name)),...
            %         'Optimization Toolbox'));
            %     if ~has_fmincon
            %         exc = SrtSetupError(['SReachSet with ''lag-over/under''',...
            %             'and bound_set_method as ''optim-box'' needs ',...
            %             'MATLAB''s Optimization Toolbox.']);
            %         throw(exc);
            %     end
            %     if any(contains(inpar.UsingDefaults,'box_center'))
            %         throw(SrtInvalidArgsError(['Box center (fixed center ',...
            %             'around which the box is optimized for) is a ',...
            %             'required input for bound_set_method: optim-box']));
            %     end
            case 'box'
                if any(contains(inpar.UsingDefaults,'err_thresh'))
                    throw(SrtInvalidArgsError(['err_thresh (threshold for ',...
                        'bisection) is a required input for ',...
                        'bound_set_method: box']));
                end
            case 'random'
                if any(contains(inpar.UsingDefaults,'num_dirs'))
                    throw(SrtInvalidArgsError(['num_dirs (no. of direction ',...
                        'vectors to create bounded polytope) is a required ',...
                        'input for bound_set_method: random']));
                end
            case 'load'
                if any(contains(inpar.UsingDefaults,'load_str'))
                    throw(SrtInvalidArgsError(['Expected a path to matfile ',...
                        'with bounded set (only) as input for ',...
                        'bound_set_method: load']));
                end
        end
    end
end



