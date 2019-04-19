function varargout = ellipsoidsFromMonteCarloSims(...
    concat_state_realization, state_dim, relv_states, plot_options)
% Get/plot (minimum volume) ellispoids corresponding to a Monte-Carlo simulation
% ============================================================================
% 
% Computes (and plots) the Lowner-John ellipsoid (Boyd and Vandenberghe, Convex
% Optimization, Ch. 8.4.1) that fits the set of samples at each time instant in
% the trajectory. This function is typically useful in conjunction with
% generateMonteCarloSims to understand the spread of the Monte-Carlo simulated
% trajectories.
%
% Usage: See examples/cwhSReachPoint.m
%
% =============================================================================
%
% set_of_ellipsoids = ellipsoidsFromMonteCarloSims(...
%    concat_state_realization, state_dim, relv_states, plot_options)
%
% Inputs:
% -------
%   concat_state_realization  
%                     - Matrix of concatenate state (row) vectors stacked
%                       columnwise. Each row comprises of the state trajectory
%                       as [x_1; x_2; ...; x_N]
%   state_dim         - State dimension
%   relv_states       - A two-element vector with indices of relevant states
%                       among the states indexed from 1 to state_dim
%   plot_options      - Plot options that are directly passed to the plot. Leave
%                       empty if plotting is not desired.
%
% Outputs:
% --------
%   set_of_ellipsoids - Set of ellipsoids that cover these points. A cell array
%                       of cells are provided with each cell containing the
%                       ellipsoid center q and its shape matrix Q. The ellipsoid
%                       is (x-q)^T Q^{-1} (x-q) <= 1.
%
% See also generateMonteCarloSims
%
% Notes:
% ------
% * Requires CVX for solving the convex optimization problem using SDPT3
% * Note that the initial state is NOT a part of the 
%   concatenated_state_realization
% * Uses code obtained from
%   http://web.cvxr.com/cvx/examples/cvxbook/Ch08_geometric_probs/min_vol_elp_finite_set.m
% 
% ============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://sreachtools.github.io/license/
% 
%

    % Parameters for the simulation
    n_mcarlo_sims = size(concat_state_realization,2);
    time_horizon = size(concat_state_realization,1)/state_dim;
    set_of_ellipsoids = cell(1,time_horizon);
    n_angles = 200;
    angles   = linspace( 0, 2 * pi, n_angles );
    
    % How many MC simulations?
    max_limit_sims = 1e2;
    if n_mcarlo_sims > max_limit_sims
        relv_mc_indx = floor(linspace(1,n_mcarlo_sims,max_limit_sims));
        relv_n_mcarlo_sims  = length(relv_mc_indx);
    else
        relv_mc_indx = 1:n_mcarlo_sims;
        relv_n_mcarlo_sims  = n_mcarlo_sims;
    end
    
    % Relevant states out of the trajectory
    relv_indx = reshape([relv_states(1):state_dim:state_dim * time_horizon;
                         relv_states(2):state_dim:state_dim * time_horizon],[],1);
    relv_concat_state_realization = concat_state_realization(relv_indx, ...
        relv_mc_indx);
    
    % Angles for the ellipsoid            
    for tindx = 1:time_horizon
        % Extract the points
        x_points = relv_concat_state_realization(2*(tindx-1)+1 : 2*tindx,:);
        
        % Problem 8.11 in CVX optimization textbook
        cvx_begin quiet
            % Gurobi can't handle this! So use SDPT3
            cvx_solver SDPT3
            variable A(2, 2) symmetric;
            variable b(2, 1);
            
            maximize( det_rootn( A ) )
            subject to
                norms(A*x_points+b*ones(1,relv_n_mcarlo_sims), 2 ) <= 1;
        cvx_end
        if ~strcmpi(cvx_status, 'Solved')
            warning('SReachTools:runtime', sprintf(['CVX failed to obtain ',...
                'the ellipsoid at %d, potentially due to numerical issues.'],...
                tindx)); 
            set_of_ellipsoids(tindx) = {{mean(x_points,2),Inf(2,2)}};
        else
            % Construct the ellipsoid by matching coeff (x-q)^T Q^{-1} (x-q) =
            % (Ax + b)^2 = (x+A\b)'A'A(x+A\b)
            Q = inv(A*A');
            q = -A\b;

            % Prepare the output
            set_of_ellipsoids(tindx) = {{q,Q}};

            if ~isempty(plot_options)
                ellipse  = A \ [ cos(angles) - b(1) ; sin(angles) - b(2) ];
                h=plot(ellipse(1,:), ellipse(2,:), plot_options{:});
                h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end
        end
    end    
    varargout{1} = set_of_ellipsoids;
end
