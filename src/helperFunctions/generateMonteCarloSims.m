function [concat_state_realization, ...
          concat_disturb_realizations, ...
          saturation_indx]= generateMonteCarloSims(...
                                               n_monte_carlo_sims, ...
                                               sys, ...
                                               initial_state, ...
                                               time_horizon, ...
                                               varargin)
% Generate Monte-Carlo simulations for a Gaussian-perturbed LTI/LTV system
% (controlled or uncontrolled)
% ============================================================================
% 
% generateMonteCarloSims produces a required number of trajectories for a
% Gaussian LTI system.
%
% See also examples/forwardStochasticReachCWH.m, examples/cwhSReachPointDemo.m
%
% =============================================================================
% [concat_state_realization, concat_disturb_realizations] = ...
%    generateMonteCarloSims(n_monte_carlo_sims, sys, initial_state, ...
%       time_horizon, optimal_input_vector, optimal_input_gain)
%
% Inputs:
% -------
%   n_monte_carlo_sims  - Number of Monte-Carlo simulation particles to be used 
%                         for estimation of the reach-avoid probability
%   sys                 - System description as a LtiSystem/LtvSystem object
%   initial_state       - Deterministic x_0
%   time_horizon        - Time horizon (N) of the stochastic reach-avoid problem
%   concat_input_vector - [Optional] Open-loop controller, a column vector of
%                         dimension (sys.input_dim*N) x 1 | Required only if 
%                         the system is controlled
%   dist_feedback_gain  - [Optional] Affine disturbance feedback gain for the
%                         concatenated disturbance vector, a matrix of dimension
%                         (sys.input_dim*N) x (sys.dist_dim*N) | Required only
%                         if the system is controlled, the controller is affine
%                         disturbance feedback, and the gain matrix must be
%                         lower block triangular (with zeros in its block
%                         diagonal elements) for causality | See Notes
%   verbose             - [Optional] Verbosity of this function when saturating
%                         affine disturbance feedback controllers
%
% Outputs:
% --------
%   concat_state_realization  - Matrix of concatenate state (row) vectors
%                               stacked columnwise. Each row comprises of the
%                               state trajectory as [x_0; x_1; x_2; ...; x_N]
%   concat_disturb_realization- Matrix of concatenate disturbance (row) vectors
%                               stacked columnwise. Each row comprises of
%                               the state trajectory as [w_0; w_1; ...; w_{N-1}]
%   saturation_indx           - Binary vector that indicates which realizations
%                               had their associated affine disturbance feedback
%                               controller saturated. Potentially non-zero only
%                               if the input_gain is non-empty
%
% Notes:
% ------
% * Assumes IID Gaussian disturbance for the LTI/LTV system. 
% * For controlled system, an open-loop controller NEEDS to be provided. The
%   optimal_input_vector should be a ((sys.input_dim) * time_horizon)-dim.
%   vector U = [u_0; u_1; ...; u_N] (column vector).
% * For uncontrolled system, the optimal_input_vector NEED NOT be provided
%   dist_feedback_gain must be lower
% * The disturbance feedback gain matrix must be lower block triangular, with
%   its block diagonal submatrices as zero matrices. This ensures that the
%   affine disturbance feedback controller's value at any point of time depends
%   only on the past disturbance values => causal controller. This function DOES
%   NOT check for this structure in the input gain. TODO
% * Affine disturbance feedback controllers CAN NOT satisfy hard control bounds
%   when the disturbance is unbounded (like Gaussian). Therefore, we will
%   saturate the controller realization (associated with the disturbance
%   realization) via projection on to the concatenated input space.
%   Specifically, we solve the corresponding optimization problem for each
%   concatenated disturbance realization W
%   
%       minimize || U - (MW + D)||_2
%       subject to 
%           U \in \mathcal{U}^T
%
%   where U is the decision variable, \mathcal{U} is the input space, T is the
%   time horizon, M is the affine disturbance feedback gain, and D is the affine
%   disturbance feedback bias.
%
% ============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% 
%

    %% Input handling 
    % Ensure that n_monte_carlo_sims is a scalar and positive
    validateattributes(n_monte_carlo_sims, {'numeric'}, {'scalar',...
        'integer', '>', 0});
    % Ensure that sys is a LtiSystem/LtvSystem object
    validateattributes(sys, {'LtiSystem','LtvSystem'}, {'nonempty'});
    % Ensure sys has a Gaussian disturbance 
    if ~strcmpi(sys.dist.type,'Gaussian')
        throwAsCaller(SrtInvalidArgsError('Expected a Gaussian-perturbed', ...
            ' LtiSystem/LtvSystem object'));
    end
    % Ensure initial state vector is valid or it is a random vector
    validateattributes(initial_state, {'RandomVector','numeric'},...
        {'nonempty'});
    switch class(initial_state)
        case 'RandomVector'
            if ~strcmpi(initial_state.type,'Gaussian') ||...
                    initial_state.Dim ~= sys.state_dim
                throwAsCaller(SrtInvalidArgsError('Expected a sys.state_dim',...
                    '-dimensional Gaussian random vector as initial_state'));
            end
        case 'numeric'
            if length(initial_state) ~= sys.state_dim
                throwAsCaller(SrtInvalidArgsError('Expected a sys.state_dim',...
                    '-dimensional initial_state'));
            end
    end
    % Ensure that time_horizon is a scalar and positive
    validateattributes(time_horizon, {'numeric'}, {'scalar',...
        'integer', '>', 0});
    
    % Input handling of optimal_input_vector and creation of
    % concat_optimal_input_vector
    concat_input_vector_repeated = zeros(time_horizon * sys.input_dim,...
        n_monte_carlo_sims);
    input_concat_dist_gain = zeros(time_horizon * sys.input_dim,...
        time_horizon * sys.dist_dim);
    saturation_indx = zeros(n_monte_carlo_sims, 1);

    verbose = 0;
    gain_present = 0;
    
    if sys.input_dim > 0 && ~isEmptySet(sys.input_space)
        if length(varargin) > 3
            throwAsCaller(SrtInvalidArgsError('Too many input arguments'));
        elseif ~isempty(varargin)
            input_concat_vector = varargin{1};
            % Ensure optimal_input_vector is a valid open loop controller
            validateattributes(input_concat_vector, {'numeric'},...
                {'nonempty','vector','size',[time_horizon*sys.input_dim 1]});
            % Repeat the same open-loop policy across all the simulations
            concat_input_vector_repeated = ...
                      repmat(input_concat_vector,1, n_monte_carlo_sims);
                  
            if length(varargin) >= 2
                gain_present = 1;
                input_concat_dist_gain = varargin{2};
                validateattributes(input_concat_dist_gain, {'numeric'},...
                    {'nonempty','size',[time_horizon * sys.input_dim,...
                             time_horizon * sys.dist_dim]});                
                if length(varargin) == 3
                    verbose = varargin{3};
                    validateattributes(verbose, {'numeric'},...
                        {'nonempty','scalar','integer','>=',0,'<=',1});
                end
            end
            
        else
            %warning('SReachTools:runTime','Setting input vectors to zero');
        end
    end

    
    % Compute the concatenated matrices for X
    [Z, H, G] = getConcatMats(sys, time_horizon);
    
    if isa(initial_state,'RandomVector')
        % Draw from the initial state PDF
        concat_initial_state = mvnrnd(initial_state.parameters.mean', ...
                                      initial_state.parameters.covariance, ...
                                      n_monte_carlo_sims)';
    else
        % Concatenation of initial state vector
        concat_initial_state = repmat(initial_state,1, n_monte_carlo_sims);
    end
    W = sys.dist.concat(time_horizon);
    % Realization of the concatenated disturbance random vectors with each
    % realization stored columnwise
    concat_disturb_realizations = mvnrnd(W.parameters.mean, ...
                                         W.parameters.covariance, ...
                                         n_monte_carlo_sims)';
    
    % Realization of the random vector X (columnwise)
    % See @LtiSystem/getConcatMats for more info on the notation used
    if verbose
        disp(['Affine disturbance feedback controller will be saturated to',...
            ' the input space via projection']);
    end

    if gain_present
        concat_input_realizations = input_concat_dist_gain *...
            concat_disturb_realizations + concat_input_vector_repeated;
        % Check how many input realizations require saturation
        [concat_input_space_A, concat_input_space_b] = getConcatInputSpace(...
            sys, time_horizon);
        concat_input_space = Polyhedron('A', concat_input_space_A,...
            'b', concat_input_space_b);
        if verbose
            fprintf(['Using Polyhedron/contains to identify realizations',...
                ' that require saturation...']);
        end
        saturation_indx = concat_input_space.contains(concat_input_realizations);        
        proj_req = find(saturation_indx == 0);
        if verbose
            disp('Done');
            fprintf('We need to saturate %d realizations.\n', length(proj_req));
        end
        
        % Saturate these realizations
        realization_counter = 0;            
        realization_frac_disp = round(length(proj_req)/5);%disp every
        for realization_indx = proj_req
            % Affine disturbance feedback controller associated with
            % concat_disturb_realization (U = MW + D)
            realization_counter = realization_counter + 1;            
            if verbose && (mod(realization_counter, realization_frac_disp)==0)
                fprintf(['Completed saturating %d/%d input realizations',...
                   ' out of %d realizations\n'], realization_counter,...
                   length(proj_req), n_monte_carlo_sims);
            end
            % Saturate the resulting inputs
            cvx_begin quiet
                variable concat_saturated_input(sys.input_dim * time_horizon, 1)
                
                minimize (norm(concat_saturated_input -...
                    concat_input_realizations(:, realization_indx)))
                subject to
                    concat_input_space_A * concat_saturated_input <=...
                        concat_input_space_b;                    
            cvx_end
            switch cvx_status
                case {'Solved','Inaccurate/Solved'}
                    % Input realization saturated via projection
                    concat_input_realizations(:, realization_indx) =...
                        concat_saturated_input;
                    if cvx_optval > 1e-8
                        saturation_indx(:, realization_indx) = 1;
                    end
                otherwise
                    throw(SrtDevError(['Saturation via projection failed. ',...
                        'This might be due to numerical issues.']));
            end
        end
    else
        concat_input_realizations = concat_input_vector_repeated;
    end
    concat_state_realization = [concat_initial_state;
        Z * concat_initial_state + H * concat_input_realizations + ...
          G * concat_disturb_realizations];
end
