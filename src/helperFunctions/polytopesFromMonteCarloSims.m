function varargout = polytopesFromMonteCarloSims(...
    concat_state_realization, state_dim, relv_states, plot_options)
% Get/plot convex hulls corresponding to a Monte-Carlo simulation
% ============================================================================
% 
% Computes (and plots) the convex hull of the Monte-Carlo simulated trajectories 
% at each time step within the horizon. This function is typically useful in 
% conjunction with generateMonteCarloSims to understand the spread of the 
% Monte-Carlo simulated trajectories.
%
% Usage: See examples/dubinsSReachPoint.m
%
% =============================================================================
%
% set_of_ellipsoids = polytopesFromMonteCarloSims(...
%    concat_state_realization, state_dim, relv_states, plot_options)
%
% Inputs:
% -------
%   concat_state_realization  
%                     - Matrix of concatenated state (column) vectors stacked 
%                       columnwise. Each column has the state trajectory 
%                       [x_1; x_2; ...; x_N]
%   state_dim         - State dimension
%   relv_states       - A two-element vector with indices of relevant states
%                       among the states indexed from 1 to state_dim
%   plot_options      - Plot options (MPT3 Polyhedron/plot) that are directly 
%                       passed to the plot. Leave empty if plotting is undesired
%
% Outputs:
% --------
%   set_of_polytopes  - Set of polytopes comprised of the points. (Array)
%   hpoly             - Handle of the MPT3 plot for one of the polytopes
%
% See also generateMonteCarloSims
%
% Notes:
% ------
% * Requires MPT3 to compute the convex hull
% * Note that the initial state is NOT a part of the 
%   concatenated_state_realization
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
    set_of_polytopes = [Polyhedron()];
    
    % How many MC simulations?
    max_limit_sims = 1e3;
    if n_mcarlo_sims > max_limit_sims
        relv_mc_indx = floor(linspace(1,n_mcarlo_sims,max_limit_sims));
    else
        relv_mc_indx = 1:n_mcarlo_sims;
    end

    % Relevant states out of the trajectory
    relv_indx = reshape([relv_states(1):state_dim:state_dim * time_horizon;
                         relv_states(2):state_dim:state_dim * time_horizon], ...
                         [],1);
    relv_concat_state_realization = concat_state_realization(relv_indx, ...
        relv_mc_indx);
    
    % Angles for the ellipsoid            
    for tindx = 1:time_horizon
        % Extract the points
        x_points = relv_concat_state_realization(2*(tindx-1)+1 : 2*tindx,:);
        
        % Prepare the output
        poly_at_tindx = Polyhedron('V', x_points');
        set_of_polytopes(tindx) = poly_at_tindx;

        if ~isempty(plot_options)
            hpoly = plot(poly_at_tindx, plot_options{:});
        end
    end    
    varargout{1} = set_of_polytopes;
    varargout{2} = hpoly;
end
