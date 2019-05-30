function varargout = SReachFwd(prob_str, sys, initial_state, target_time, ...
    varargin)
% Perform forward stochastic reachability analysis of a linear system with
% additive stochastic disturbance, initialized to a random vector
% ============================================================================
%
% Forward stochastic reachability analysis quantifies the stochasticity of the
% state at a future time instant, as well as the probability of the state lying
% in a pre-specified set at a future time instant. This function implements
% ideas from
%
% A. Vinod, B. HomChaudhuri, and M. Oishi, "Forward Stochastic Reachability
% Analysis for Uncontrolled Linear Systems using Fourier Transforms", In
% Proceedings of the 20th International Conference on Hybrid Systems:
% Computation and Control (HSCC), 2017.
%
% While the Fourier transform-based results apply for arbitrary 
% distributions, SReachFwd considers only Gaussian-perturbed LTI systems. 
% In this case, the approach coincides with Kalman filter updates, and it 
% is grid-free, recursion-free, and sampling-free. 
%
% SReachFwd also exploits the functionality of random vector to provide
% forward stochastic reachability analysis of linear systems perturbed by
% non-Gaussian disturbance. In this case, we use Monte-Carlo simulations to
% estimate the mean, covariance, and the probability of lying in a set.
%
% See also examples/cwhSReachFwd.m.
%
% ============================================================================
% 
% varargout = SReachFwd(prob_str, sys, initial_state, target_time, varargin)
% 
% Inputs:
% -------
%   prob_str      - String specifying the problem of interest
%                       1. 'state-stoch' : Provide mean and covariance at state
%                                          at the spjecified time
%                       2. 'state-prob'  : Compute the probability that the
%                                          state will lie in a polytope at the
%                                          specified time
%                       3. 'concat-stoch': Provide mean and covariance of the
%                                          concatenated state vector up to a
%                                          specified time
%                       4. 'concat-prob' : Compute the probability that the
%                                          concatenated state vector up to a
%                                          specified time lies in the given
%                                          target tube
%   sys           - System description as a LtiSystem/LtvSystem object
%   initial_state - Initial state as a deterministic n-dimensional vector
%                   or a RandomVector object
%   target_time   - Time of interest (non-negative scalar) | If target_time = 0,
%                   then the stochasticity of the initial state is analyzed.
%   target_set/tube  
%                 - [Required only for state/concat-prob] Polyhedron/Tube object
%                   over which the probability must be computed
%   desired_accuracy 
%                 - [Optional for state/concat-prob] Maximum absolute deviation
%                   from the true probability estimate [Default: 1e-2]
%
% Outputs:
% --------
%   rv            - ['state/concat-stoch'] Random vector describing the
%                   state / concatenated state vector | It is a vector of
%                   dimension (target_time + 1) * sys.state_dim to include the
%                   initial state
%   prob          - ['state/concat-prob'] Probability of occurence
%
% Notes:
% ------
% * Assumes IID disturbance.
% * The outputs are either (mean_vec, cov_mat) or (prob), depending on the
%   method_str
% * For concat-prob, target_time can be any value between 0 and N where
%   target_tube has a length of N+1 (N+1 target sets to include the
%   constraints on the initial state).
% * For XXX-prob, the random vector is provided as the second output.
% * Prob{ w \in \theta Polytope(A,b) } is computed using
%   RandomVector/getProbPolyhedron. Note that this is an underapproximation of
%   the true probability
% * Requires the desired_accuracy to be at least 1e-2.
%
% ============================================================================
%
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://sreachtools.github.io/license/
%
%

    % Input parsing
    valid_prob_str = {'state-stoch','state-prob','concat-stoch','concat-prob'};
    inpar = inputParser();
    inpar.addRequired('prob_str', @(x) any(validatestring(x,valid_prob_str)));
    inpar.addRequired('sys', @(x) validateattributes(x, ...
        {'LtiSystem','LtvSystem'}, {'nonempty'}));
    inpar.addRequired('initial_state', @(x) validateattributes(x, ...
        {'RandomVector', 'numeric'}, {'nonempty'}))
    inpar.addRequired('target_time', @(x) validateattributes(x, ...
        {'numeric'}, {'scalar', 'integer', '>=', 0}));

    try
        inpar.parse(prob_str, sys, initial_state, target_time);
    catch err
        exc = SrtInvalidArgsError.withFunctionName();
        exc = exc.addCause(err);
        throwAsCaller(exc);
    end

    % Decide the approach to take
    prob_str_splits = split(prob_str, '-');
    
    % Ensure that:
    % 1. Initial state is a column vector of dimension sys.state_dim OR
    %    a RandomVector (Gaussian) object of dimension sys.state_dim
    % 2. Given system is an uncontroller LTI/LTV system with Gaussian 
    %    disturbance
    % 3. For prob computation, ensure the optional arguments are all ok
    otherInputHandling(sys, initial_state, prob_str_splits, varargin, ...
        target_time);

    if target_time > 0
        % IID assumption allows to compute the mean and covariance of the
        % concatenated disturbance vector W
        concat_disturb = sys.dist.concat(target_time);

        % Compute the state_trans_matrix and controllability matrix for
        % disturbance | No input
        [Z,~,G] = sys.getConcatMats(target_time);

        if strcmpi(prob_str_splits{1},'state')
            state_trans_mat = Z(end-sys.state_dim+1:end, ...
                                        end-sys.state_dim+1:end);
            flipped_ctrb_mat_disturb = G(end-sys.state_dim+1:end,:);

            rv = state_trans_mat * initial_state + flipped_ctrb_mat_disturb * ...
                concat_disturb;
        elseif strcmpi(prob_str_splits{1},'concat')
            rv = Z * initial_state + G * concat_disturb;
            if isa(initial_state, 'RandomVector')
                rv = [initial_state;
                      rv];
            else
                rv = [initial_state;zeros(rv.dim,1)] + ...
                    [zeros(sys.state_dim, rv.dim);
                     eye(rv.dim)] * rv;                
            end
        end
    elseif isa(initial_state, 'RandomVector')
        rv = initial_state;
    else
        throw(SrtInvalidArgsError('Initial state is not random'));
    end
    
    if strcmpi(prob_str_splits{2},'prob')
        if length(varargin) == 2
            desired_accuracy = varargin{2};
            validateattributes(desired_accuracy, {'numeric'}, {'scalar', ...
                '>=', 1e-2, '<', 1});
        elseif length(varargin) == 1
            desired_accuracy = 1e-2;
        else
            throwAsCaller(SrtInvalidArgsError('Too many input arguments'));
        end
        if isa(rv, 'numeric')
            target_set = varargin{1};                
            % Not really a random vector
            prob = target_set.contains(rv);
        elseif strcmpi(prob_str_splits{1},'state')
            target_set = varargin{1};                      
            % Compute probability at time target_time of x \in target_set
            prob = rv.getProbPolyhedron(target_set, desired_accuracy);
        elseif strcmpi(prob_str_splits{1},'concat')
            target_tube = varargin{1};
            [concat_safety_tube_A, concat_safety_tube_b] = ...
                target_tube.concat([1 target_time+1]);
            polytope_for_concat_tube = Polyhedron('H', ...
                [concat_safety_tube_A, concat_safety_tube_b]);
            prob = rv.getProbPolyhedron(polytope_for_concat_tube, ...
                desired_accuracy);
        end
        varargout{1} = prob;
        varargout{2} = rv;
    elseif strcmpi(prob_str_splits{2},'stoch')
        varargout{1} = rv;
    end
end

function otherInputHandling(sys, initial_state, prob_str_splits, ...
    optional_args, target_time)
    % Ensure that initial state is a column vector of appropriate dimension OR
    % a Gaussian random vector of approriate dimension 
    if isa(initial_state,'RandomVector') 
        if strcmpi(initial_state.type, 'Gaussian') &&...
            initial_state.dim==sys.state_dim
        else
            throwAsCaller(SrtInvalidArgsError(['Expected a sys.state_dim-',...
                'dimensional Gaussian random vector for initial state']));
        end
    elseif isa(initial_state,'numeric') &&...
            ~isequal(size(initial_state), [sys.state_dim 1])
        throwAsCaller(SrtInvalidArgsError(['Expected a sys.state_dim-', ...
            'dimensional column-vector for initial state']));
    end

    % Ensure that the given system has a Gaussian disturbance
    if isa(sys.dist, 'RandomVector') && ~strcmpi(sys.dist.type, 'Gaussian')
        throwAsCaller(SrtInvalidArgsError('Expected a Gaussian-perturbed ', ...
            'LTI/LTV system'));
    elseif ~isa(sys.dist, 'RandomVector')
        throwAsCaller(SrtInvalidArgsError('Expected a stochastic system'));
    end

    % Ensure that the given system is uncontrolled
    if sys.input_dim ~= 0
        throwAsCaller(SrtInvalidArgsError('Expected an uncontrolled system'));
    end
    
    % Ensure the optional arguments for prob computation are all ok
    if strcmpi(prob_str_splits{2},'prob')
        if length(optional_args)~=2
            throwAsCaller(SrtInvalidArgsError(['Expected {target set/target',...
                ' tube} and desired accuracy']));
        end   
        % Ensure target_set is a non-empty Polyhedron
        switch prob_str_splits{1}
            case 'state'
                target_set = optional_args{1};
                % Ensure target_set is a non-empty Polyhedron
                if ~(isa(target_set, 'Polyhedron') && ...
                     ~target_set.isEmptySet() && ...
                     target_set.Dim == sys.state_dim)

                    throwAsCaller(SrtInvalidArgsError(['Expected a non-', ...
                        'empty polyhedron of dimension sys.state_dim as ',...
                        'target set']));
                end
            case 'concat'
                target_tube = optional_args{1};
                if ~(isa(target_tube, 'Tube') && ...
                     target_tube.tube(1).Dim == sys.state_dim && ...
                     length(target_tube) >= target_time + 1) 

                    throwAsCaller(SrtInvalidArgsError(['Expected a target ', ...
                        'tube of length not smaller than target_time+1 and ',...
                        'dimension sys.state_dim']));
                end
        end
    end    
end
