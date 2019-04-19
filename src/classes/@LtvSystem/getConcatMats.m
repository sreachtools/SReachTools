function [Z,H,G] = getConcatMats(sys, time_horizon)
% Get concatenated matrices
% ============================================================================
% 
% Computes the matrices corresponding to the concatentated state vector X.
%
% Consider a LtvSystem object with n as the state_dim, m as the
% input_dim, and p as the disturbance_dim. Given a time of
% interest N, we define a concatenated state vector (a nN-dimensional vector)
%           __       __
%           |   x_1   |
%           |   x_2   |
%       X = |   ...   |
%           | x_{N-1} |          
%           |  x_{N}  |          
%           ---     ---
% where x_t is the state of the system with 1 <= t <= N.  Similarly, one can
% define concated input and noise vectors U and W (mN-dimensional and
% pN-dimensional vectors),
%           __       __         __       __
%           |   u_0   |         |   w_0   |
%           |   u_1   |         |   w_1   |
%       U = |   ...   |,   W  = |   ...   |
%           | u_{N-2} |         | w_{N-2} |      
%           | u_{N-1} |         | w_{N-1} |      
%           ---     ---         ---     ---
%
% Given the initial state x_0, we have
%
%       X = Z * x_0 + H * U + G * W
%
% where Z (nN x n matrix), H (nN x mN matrix), and G  (nN x
% pN matrix) are appropriate matrices. These matrices (with minor
% modifications noted below) are given in (3) in 
%    J. Skaf and S. Boyd, "Design of Affine Controllers via Convex
%    Optimization", in IEEE Trans. Automatic Control, 2010. 
%
% This function computes Z, H, and G.
%
% Usage:
% ------
%
% % Compute the concatenated matrices for a double integrator with a time of
% % interest, 10
%
% % Problem parameters
% time_horizon = 10;
% T = 0.25;
% umax = 0.75;
% dmax = 0.1;
% % Double integrator system
% sys = LtvSystem(...
%     'StateMatrix', [1, T; 0, 1], ...
%     'InputMatrix', [T^2; T], ...
%     'InputSpace', Polyhedron('lb', -umax, 'ub', umax), ...
%     'DisturbanceMatrix', eye(2), ...
%     'Disturbance', Polyhedron('lb', -dmax *ones(2,1), 'ub', dmax *ones(2,1)));
% 
% % Get the concatenated matrices
% [Z,H,G] = getConcatMats(sys, time_horizon);
%
% =============================================================================
%
% [Z,H,G] = getConcatMats(sys, time_horizon)
%
% Inputs:
% -------
%   sys          - An object of LtvSystem class 
%   time_horizon - Time horizon (N) with the control provided from 0 to N-1
%
% Outputs:
% --------
%   Z - Concatenated state matrix
%   H - Concatenated input matrix
%   G - Concatenated disturbance matrix
%
% Notes:
% ------
% * For control-free and/or disturbance-free LTI/LTV systems, H and G are set to
%   zeros( sys.state_dim * time_horizon, 1) as appropriate.
% * Deviation from Skaf and Boyd's definition,
%     * Concatenated state is X=[x_1 x_2 ... x_{N}], with the initial state
%     x_0 EXCLUDED in contrast to Skaf and Boyd's TAC 2010 definitions.
% * Computes the extended controllability matrix via for loops. (suboptimal way)
% * This function also serves as a delegatee for input handling
% 
% ============================================================================
%
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://sreachtools.github.io/license/
% 
%

    % Ensure that time_horizon is a scalar
    assert( isscalar(time_horizon) && time_horizon > 0, ...
           'SReachTools:invalidArgs', ...
           'Expected a scalar positive time_horizon');

    %% Construct Z matrix --- concatenated state matrix
    % Z = [A;A^2;...A^{N}]
    Z = zeros(sys.state_dim * time_horizon, sys.state_dim);
    for t_indx = 1:time_horizon
        Z((t_indx-1)*sys.state_dim +1 : t_indx*sys.state_dim,:) =...
            computeMatProductsStateMat(sys, 0, t_indx);
    end

    %% Construct H matrix --- concatenated input matrix
    if sys.input_dim > 0
        H = zeros(sys.state_dim * time_horizon, sys.input_dim * time_horizon);
        for t_indx = 1:time_horizon
            row_for_H = zeros(sys.state_dim, sys.input_dim * time_horizon);
            for sub_indx = 1:t_indx
                if sys.islti()
                    input_matrix_temp = sys.input_mat;
                else
                    input_matrix_temp = sys.input_mat(sub_indx-1);
                end
                row_for_H(:, ...
                    (sub_indx-1)*sys.input_dim + 1: sub_indx*sys.input_dim) =...
                        computeMatProductsStateMat(sys, sub_indx, t_indx) * ...
                        input_matrix_temp;
            end            
            H((t_indx-1)*sys.state_dim +1 : t_indx*sys.state_dim,:) = row_for_H;
        end
    else
        H = zeros(sys.state_dim * time_horizon, 0);
    end

    %% Construct G matrix --- concatenated disturbance matrix
    if sys.dist_dim > 0
        G = zeros(sys.state_dim * time_horizon, sys.dist_dim * time_horizon);
        for t_indx = 1:time_horizon
            row_for_G = zeros(sys.state_dim, sys.dist_dim * time_horizon);
            for sub_indx = 1:t_indx
                if sys.islti()
                    dist_matrix_temp = sys.dist_mat;
                else
                    dist_matrix_temp = sys.dist_mat(sub_indx-1);
                end
                row_for_G(:, ...
                    (sub_indx-1)*sys.dist_dim + 1: sub_indx*sys.dist_dim) =...
                        computeMatProductsStateMat(sys, sub_indx, t_indx) * ...
                        dist_matrix_temp;
            end            
            G((t_indx-1)*sys.state_dim +1 : t_indx*sys.state_dim,:) = row_for_G;
        end
    else
        G = zeros(sys.state_dim * time_horizon, 0);
    end
end

function prod_state_mat = computeMatProductsStateMat(sys, t, tau)
    % To handle time-varying setup
    prod_state_mat = eye(sys.state_dim);
    if sys.islti()
        if t<tau
            prod_state_mat = sys.state_mat^(tau-t);
        else
            % Return identity
        end
    else
        if t<tau
            for t_indx = t:tau-1
                prod_state_mat = prod_state_mat * sys.state_mat(t_indx);
            end    
        else
            % Return identity
        end
    end
end
