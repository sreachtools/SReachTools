function varargout = SReachFwd(prob_str, sys, initial_state, target_time, ...
    varargin)
% Perform forward stochastic reachability analysis of a Gaussian-perturbed
% linear system
% ============================================================================
%
% Perform forward stochastic reachability analysis of a Gaussian-perturbed
% linear system. This function implements ideas from
%
% A. Vinod, B. HomChaudhuri, and M. Oishi, "Forward Stochastic Reachability
% Analysis for Uncontrolled Linear Systems using Fourier Transforms", In
% Proceedings of the 20th International Conference on Hybrid Systems:
% Computation and Control (HSCC), 2017.
%
% Usage: See examples/forwardStochasticReachCWH.mlx.
%
% ============================================================================
% 
% varargout = SReachFwd(prob_str, sys, initial_state, target_time, varargin)
% 
% Inputs:
% -------
%   prob_str      - String specifying the problem of interest
%                       1. 'state-stoch' : Provide mean and covariance at state
%                                          at the specified time
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
%   target_time   - Time of interest (positive scalar)
%   target_set/tube 
%                 - [Required only for state/concat-prob] Polytope/Tube
%                   over which the probability must be computed
%   desired_accuracy
%                 - [Required only for state/concat-prob] Accuracy for the
%                   integral
%   
%
% Outputs:
% --------
%   mean_vec      - ['state/concat-stoch'] Mean of the stochastic disturbance
%   cov_mat       - ['state/concat-stoch'] Covariance of the stochastic disturbance
%   prob          - ['state/concat-prob'] Probability of occurence
%
% Notes:
% ------
% * Requires Gaussian-perturbed LTI/LTV system
% * Assumes IID disturbance.
% * Few minor components are not covered by the unit tests (See TODO-Test)
%
% ============================================================================
%
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
%
%

    % Input parsing
    valid_prob_str = {'state-stoch','state-prob','concat-stoch','concat-prob'};
    inpar = inputParser();
    inpar.addRequired('prob_str', @(x) any(validatestring(x,valid_prob_str)));
    inpar.addRequired('sys', @(x) validateattributes(x,...
        {'LtiSystem','LtvSystem'}, {'nonempty'}));
    inpar.addRequired('initial_state', @(x) validateattributes(x, ...
        {'RandomVector', 'numeric'}, {'nonempty'}))
    inpar.addRequired('target_time', @(x) validateattributes(x, ...
        {'numeric'}, {'scalar', 'integer', '>', 0}));

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
    % 2. Given system is an uncontroller LTI/LTV system with Gaussian disturbance
    % 3. For prob computation, ensure the optional arguments are all ok
    otherInputHandling(sys, initial_state, prob_str_splits, varargin, target_time);

    % IID assumption allows to compute the mean and covariance of the
    % concatenated disturbance vector W
    mean_concat_disturb = kron(ones(target_time,1), sys.dist.parameters.mean);
    cov_concat_disturb =kron(eye(target_time), sys.dist.parameters.covariance);

    % Compute the state_trans_matrix and controllability matrix for
    % disturbance
    [Z,~,G] = sys.getConcatMats(target_time);

    if strcmpi(prob_str_splits{1},'state')
        state_trans_mat = Z(end-sys.state_dim+1:end,...
                                    end-sys.state_dim+1:end);
        flipped_ctrb_mat_disturb = G(end-sys.state_dim+1:end,:);

        if isa(initial_state,'RandomVector')
            % mu_x = A^\tau * mu_{x_0} + C * mu_{W};
            mean_vec = state_trans_mat * initial_state.parameters.mean + ...
                           flipped_ctrb_mat_disturb * mean_concat_disturb;
            % cov_x = A^\tau * cov_{x_0} * (A^\tau)' + C * cov_{W} * C';
            cov_mat =state_trans_mat* initial_state.parameters.covariance *...
                    state_trans_mat' + flipped_ctrb_mat_disturb *...
                                     cov_concat_disturb * flipped_ctrb_mat_disturb';
        else
            % mu_x = A^\tau * x_0 + C * mu_W;
            mean_vec = state_trans_mat * initial_state + ...
                                    flipped_ctrb_mat_disturb * mean_concat_disturb;
            % cov_x = C * cov_{W} * C';
            cov_mat = flipped_ctrb_mat_disturb * cov_concat_disturb *...
                                                          flipped_ctrb_mat_disturb';
        end
    elseif strcmpi(prob_str_splits{1},'concat')
        if isa(initial_state,'RandomVector')
            % mu_X = Z * mu_{x_0} + G * mu_{W};
            mean_vec = Z * initial_state.parameters.mean + G * mean_concat_disturb;
            % cov_X = Z * cov_{x_0} * Z' + G * cov_{W} * G';
            cov_mat = Z * initial_state.parameters.covariance * Z' +...
                        G * cov_concat_disturb * G';
        else
            % mu_X = Z * x_0 + G * mu_{W};
            mean_vec = Z * initial_state + G * mean_concat_disturb;
            % cov_X = G * cov_{W} * G';
            cov_mat = G * cov_concat_disturb * G';
        end
    end
    if ~issymmetric(cov_mat)
        % Compute the symmetric component of it
        symm_cov_mat = (cov_mat+cov_mat')/2;
        % Max error element-wise
        max_err = max(max(abs(cov_mat - symm_cov_mat)));
        if max_err > eps
            % TODO-Test: Not sure when the matrices are non-symmetric
            warning('SReachTools:runtime',sprintf(['Non-symmetric ',...
                'covariance matrix made symmetric (max element-wise error:',...
                '%1.3e)!'], max_err));
        end
        cov_mat = symm_cov_mat;
    end
    if strcmpi(prob_str_splits{2},'prob')
        desired_accuracy = varargin{2};
        
        if strcmpi(prob_str_splits{1},'state')
            % Compute probability at time target_time of x \in target_set
            target_set = varargin{1};

            % Construct the half-space representation for qscmvnv
            qscmvnv_lb = repmat(-Inf, [size(target_set.A, 1), 1]);
            qscmvnv_coeff_matrix = target_set.A;
            qscmvnv_ub = target_set.b - target_set.A * mean_vec;

        elseif strcmpi(prob_str_splits{1},'concat')
            target_tube = varargin{1};

            % Get half space representation of the target tube until the
            % target_time of interest (guaranteed to be smaller than the length
            % of the target tube)
            [concat_target_tube_A, concat_target_tube_b] =...
                target_tube.concat([2 target_time+1]);

            % Construct the half-space representation for qscmvnv
            qscmvnv_lb  = repmat(-Inf, [size(concat_target_tube_A, 1), 1]);
            qscmvnv_coeff_matrix = concat_target_tube_A;
            qscmvnv_ub  = concat_target_tube_b -concat_target_tube_A * mean_vec;

        end
        % Call Genz's algorithm in an iterative approach to compute the
        % probability. Uses the desired_accuracy and the error_estimate from
        % qscmvnv to navigate the number of particles used in qscmvnv
        try
            prob = iteratedQscmvnv(cov_mat, ...
                                   qscmvnv_lb, ...
                                   qscmvnv_coeff_matrix, ...
                                   qscmvnv_ub, ...
                                   desired_accuracy, ...
                                   10);
        catch
            %TODO-Test: Not sure when qscmvnv will bug out!
            err=SrtDevError(['Error in qscmvnv (Quadrature of ',...
                'multivariate Gaussian)']);
            throw(err);
        end
        varargout{1} = prob;
    elseif strcmpi(prob_str_splits{2},'stoch')
        varargout{1} = mean_vec;
        varargout{2} = cov_mat;
    end
end

function otherInputHandling(sys, initial_state, prob_str_splits, optional_args, target_time)
    % Ensure that initial state is a column vector of appropriate dimension OR
    % random vector of approriate dimension
    if isa(initial_state,'RandomVector') && initial_state.dim~=sys.state_dim...
            && strcmp(initial_state.type, 'Gaussian')
        %TODO-Test: no support for non-Gaussian as of yet
        err = SrtInvalidArgsError(['Expected a sys.state_dim-dimensional ',...
            'Gaussian random vector for initial state']);
        throw(err);
    elseif isa(initial_state,'numeric') &&...
            ~isequal(size(initial_state), [sys.state_dim 1])
        err = SrtInvalidArgsError(['Expected a sys.state_dim-dimensional ',...
            'column-vector for initial state']);
        throw(err);
    end

    % Ensure that the given system has a Gaussian disturbance
    if isa(sys.dist, 'RandomVector') && ~strcmpi(sys.dist.type, 'Gaussian')
        %TODO-Test: no support for non-Gaussian as of yet
        err=SrtInvalidArgsError('Expected a Gaussian-perturbed LTI/LTV system');
        throw(err);
    elseif ~isa(sys.dist, 'RandomVector')
        err = SrtInvalidArgsError('Expected a stochastic LTI/LTV system');
        throw(err);
    end

    % Ensure that the given system is uncontrolled
    if sys.input_dim ~= 0
        err = SrtInvalidArgsError('Expected an uncontrolled LTI/LTV systems');
        throw(err);
    end
    
    % Ensure the optional arguments for prob computation are all ok
    if strcmpi(prob_str_splits{2},'prob')
        if length(optional_args)~=2
            err=SrtInvalidArgsError(['Expected a target set or target tube ',...
                'and a desired accuracy']);
            throw(err);
        end   
        desired_accuracy = optional_args{2};
        validateattributes(desired_accuracy, {'numeric'}, {'scalar', '>', 0})
        % TODO-Test: Triggers costly computations in testing
        if sys.state_dim <=3 && desired_accuracy < 1e-8
            warning('SReachTools:desiredAccuracy',...
                'desired_accuracy < 1e-8 might be hard to enforce');
        elseif sys.state_dim >3 && desired_accuracy < 1e-4
            warning('SReachTools:desiredAccuracy',...
                'desired_accuracy < 1e-4 might be hard to enforce');
        end
        % Ensure target_set is a non-empty Polyhedron
        switch prob_str_splits{1}
            case 'state'
                target_set = optional_args{1};
                % Ensure target_set is a non-empty Polyhedron
                if ~(isa(target_set, 'Polyhedron') && ~target_set.isEmptySet()... 
                    && target_set.Dim == sys.state_dim)
                    err=SrtInvalidArgsError(['Expected a non-empty polyhedron ',...
                        'of dimension sys.state_dim as target set']);
                    throw(err);
                end
            case 'concat'
                target_tube = optional_args{1};
                if ~(isa(target_tube, 'Tube') &&...
                    target_tube.tube(1).Dim == sys.state_dim &&...
                    length(target_tube) >= target_time + 1) 
                    err=SrtInvalidArgsError(['Expected a target tube of length',...
                        ' not smaller than target_time+1 and dimension ',...
                        ' sys.state_dim']);
                    throw(err);
                end     
        end
    end    
end
