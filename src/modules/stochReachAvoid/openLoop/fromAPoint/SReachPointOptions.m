function options = SReachPointOptions(prob_str, method_str)
    options.prob_str = prob_str;
    options.method_str = method_str;            
    v = ver;            
    switch lower(method_str)
        case 'genzps-open'            
            options.desired_accuracy = 1e-3;
            options.PSoptions = psoptimset('display','off');
            has_patternsearch = any(strcmp(cellstr(char(v.Name)),...
                'Global Optimization Toolbox'));
            if ~has_patternsearch
                exc = SrtSetupError(['SReachPoint with ''genzps-open'' ',...
                    'option needs MATLAB''s Global Optimization Toolbox.']);
                throw(exc);
            end
        case 'chance-open'
            options.desired_accuracy = 1e-3;                
        case 'chance-affine'
            options.desired_accuracy = 1e-2;
            options.max_input_viol_prob = 1e-2;
            options.bisect_lb = 0;
            options.bisect_ub = 1;
            options.tau_initial = 1;     % Weight for the DC constraint viol
            options.scaling_tau = 2;     % Multiplying factor for the weight
            options.tau_max = 1e5;       % Saturation for the weight to avoid numerical issues
            options.iter_max = 20;       % Max number of DC iterations
            options.dc_conv_tol = 1e-4;  % DC convergence threshold
            options.slack_tol = 1e-8;    % When is the slack variable == to the norm values?
            options.verbose = 1;
    end
end

%   guess_optimal_input_vector
%                    - (Optional) Provide a concatenated guess for the optimal
%                      input policy vector in the form of U = [u_0; u_1; ...;
%                      u_N] for patternsearch. [Default []. This will trigger a
%                      CVX-based computation for initialization of
%                      patternsearch.]
%   desired_accuracy - (Optional) Accuracy for patternsearch  [Default 5e-3]
%   PSoptions        - (Optional) Options for patternsearch [Default
%                       psoptimset('Display', 'off')]

% * Method 'genzps' has the following dependencies
%       * EXTERNAL DEPENDENCY: Uses CVX (optional)
%                              Needs CVX to setup a convex optimization problem
%                              that initializes the patternsearch-based
%                              optimization. If CVX is unavailable, the user may
%                              provide a guess for the initialization.
%       * Specify both desired_accuracy and PSoptions or neither to use the
%         defaults 
%       * Specify an optional guess_optimal_input_vector to skip the use of CVX
% * Method 'cccpwl' has the following dependencies
%       * EXTERNAL DEPENDENCY: Uses CVX 
%                              Needs CVX to setup the convex chance-constrained
%                              problems
%       * Specify a desired_accuracy if required. Else, a default value of 1e-3
%           is used.
% * Method 'ccciter' has the following dependencies
%       * EXTERNAL DEPENDENCY: Uses CVX 
%                              Needs CVX to setup the convex chance-constrained
%                              problems
%       * Specify a desired_accuracy if required. Else, a default value of 1e-3
%           is used.

        % Stochastic reach-avoid probability may be non-trivial
        % EXTERNAL DEPENDENCY CHECK
% 
% 
% 
% 
%             case 'genzps-open'
%                 % Parsing the optional arguments 
%                 if length(varargin) == 3
%                     % First optional argument is the guess_opt_controller
%                     if ~isempty(varargin{1})
%                         assert(isvector(varargin{1}) &&...
%                                length(varargin{1}) == ...
%                                sys.input_dim * time_horizon, ...
%                                'SReachTools:invalidArgs', ...
%                                ['Expected a well-dimensioned row vector ', ...
%                                 'guess_opt_controller']);
%                         guess_opt_controller = varargin{1};
%                     else
%                         guess_opt_controller = [];
%                     end
%                     % Second optional argument is the desired_accuracy
%                     assert(isscalar(varargin{2}), ...
%                            'SReachTools:invalidArgs', ...
%                            'Expected a scalar value for desired_accuracy');
%                     desired_accuracy = varargin{2};
%                     % Third optional argument is the options for patternsearch,
%                     % PSoptions (TODO: No validation being done here)
%                     PSoptions = varargin{3};
%                 elseif length(varargin) == 1
%                     % First optional argument is the guess_opt_controller
%                     if ~isempty(varargin{1})
%                         assert(isvector(varargin{1}) &&...
%                                length(varargin{1}) == ...
%                                sys.input_dim * time_horizon, ...
%                                'SReachTools:invalidArgs', ...
%                                ['Expected a well-dimensioned row vector ', ...
%                                 'guess_opt_controller']);
%                         guess_opt_controller = varargin{1};
%                     else
%                         guess_opt_controller = [];
%                     end
%                     desired_accuracy = 1e-3;
%                     PSoptions = psoptimset('Display', 'off');
%                 elseif isempty(varargin)
%                     guess_opt_controller = [];
%                     desired_accuracy = 1e-3;
%                     PSoptions = psoptimset('Display', 'off');
%                 else
%                     exc = SrtInvalidArgsError([...
%                         'guess_opt_controller, desired_accuracy, ', ...
%                         'and PSoptions are the only additional options.']);
%                     throw(exc);
%                 end
% 
%             case 'chance-open'
%                 if length(varargin) == 1
%                     % Optional argument is the desired_accuracy
%                     assert(isscalar(varargin{1}), ...
%                            'SReachTools:invalidArgs', ...
%                            'Expected a scalar value for desired_accuracy');
%                     desired_accuracy = varargin{1};
%                 elseif isempty(varargin)
%                     desired_accuracy = 1e-3;
%                 else
%                     exc = SrtInvalidArgsError(['More inputs provided to ',...
%                            'getLowerBoundStochasticReachAvoid than expected.']);
%                     throw(exc);
%                 end
%             case 'chance-affine'
%                 switch length(varargin)
%                     case 3
%                         desired_accuracy = varargin{1};
%                         max_state_violation_prob = varargin{2};
%                         max_input_violation_prob = varargin{3};
%                     case 2
%                         desired_accuracy = varargin{1};
%                         max_state_violation_prob = varargin{2};
%                         max_input_violation_prob = 0.001;
%                     case 1
%                         desired_accuracy = varargin{1};
%                         max_state_violation_prob = 0.1;
%                         max_input_violation_prob = 0.001;
%                     case 0
%                         desired_accuracy = 5e-3;
%                         max_state_violation_prob = 0.1;
%                         max_input_violation_prob = 0.001;
%                     otherwise
%                         exc = SrtInvalidArgsError(['More inputs provided to ',...
%                                'getLowerBoundStochasticReachAvoid than expected.']);
%                         throw(exc);
%                 end
%         % EXTERNAL DEPENDENCY CHECK --- CVX
%         assert(exist('cvx_begin','file')==2, ...
%                'SReachTools:setup_error', ...
%                ['This function uses CVX if no initial guess (input ', ...
%                 'policy) for the solver is provided. Please get CVX from ', ...
%                 'http://cvxr.com.']);
