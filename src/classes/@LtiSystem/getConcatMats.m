function [Z,H,G] = getConcatMats(sys, time_horizon)
% SReachTools/LtiSystem/getConcatMats: Get concatenated matrices
% ============================================================================
% 
% Computes the matrices corresponding to the concatentated state vector X.
%
% Consider a LtiSystem object with n as the state_dim, m as the
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
% sys = LtiSystem(...
%     'StateMatrix', [1, T; 0, 1], ...
%     'InputMatrix', [T^2; T], ...
%     'InputSpace', Polyhedron('lb', -umax, 'ub', umax), ...
%     'DisturbanceMatrix', eye(2), ...
%     'Disturbance', Polyhedron('lb', -dmax *ones(2,1), 'ub', dmax *ones(2,1)));
% % Compute the robust reach-avoid set
% [Z,H,G] = getConcatMats(sys, time_horizon);
%
% =============================================================================
%
% [Z,H,G] = getConcatMats(sys, time_horizon)
% Inputs:
% -------
%   sys          - An object of LtiSystem class 
%   time_horizon - Time of interest (N)
%
% Outputs:
% --------
%   Z - Concatenated state matrix
%   H - Concatenated input matrix
%   G - Concatenated disturbance matrix
%
% Notes:
% ------
% * For control-free and/or disturbance-free LTI systems, H and G are set to
%   zeros( sys.state_dim * time_horizon, 1) as appropriate.
% * Deviation from Skaf and Boyd's definition,
%     * Concatenated state is X=[x_1 x_2 ... x_{N}].
%     * Z definition excludes the initial state x_0 in contrast to Skaf and
%       Boyd's definition of x_0.
%     * H, G does include the first row since initial state is not there.
%     * This function computes for a LTI system instead of the original LTV
%       formulation.
% * Computes the extended controllability matrix via for loops. (suboptimal way)
% * This function also serves as a delegatee for input handling
% 
% ============================================================================
%
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/abyvinod/SReachTools/blob/master/LICENSE
% 
%

    % Ensure that time_horizon is a scalar
    assert( isscalar(time_horizon) && time_horizon > 0, ...
           'SReachTools:invalidArgs', ...
           'Expected a scalar positive time_horizon');

    %% Construct Z matrix --- concatenated state matrix
    % Z = [A;A^2;...A^{N}]
    Z=[];
    for time_index=1:time_horizon
        Z=[Z;
           sys.state_mat^time_index];
    end
    
    %% Construct H matrix --- concatenated input matrix
    if sys.input_dim > 0
        % Compute the respective extended controllability matrix (flipped)
        % [A^{N-1}B A^{N-2}B ... AB B]
        flipped_controllability_matrix_input = sys.input_mat;
        for time_index=1:time_horizon-1
            % Prepend A times (A^(time_index-1) * B) to the existing 
            % flipped_controllability_matrix_input
            flipped_controllability_matrix_input = ...
                [sys.state_mat * ...
                    flipped_controllability_matrix_input(:, ...
                                                    1:sys.input_dim), ...
                 flipped_controllability_matrix_input];
        end
        % use parts of flipped_controllability_matrix_input to obtain H 
        % blocks_of_zero_required are the number of columns in the zero matrices
        %   that have state_dim number of rows. These zero matrices go to
        %   the right of the submatrices of flipped_controllability_matrix_input
        %   in the H construction.
        % Start H construction from the top row
        H = [];
        for blocks_of_zero_required = time_horizon-1:-1:0
            relevant_indices_for_the_current_row = ...
                      ((blocks_of_zero_required * sys.input_dim) + 1) :... 
                      (sys.input_dim * time_horizon);
            current_row = [ flipped_controllability_matrix_input(:, ...
                                relevant_indices_for_the_current_row), ...
                            zeros(sys.state_dim, ...
                                blocks_of_zero_required * sys.input_dim)];
            H = [H;
                 current_row];
        end
    else
        H = zeros(sys.state_dim * time_horizon, 1);
    end

    %% Construct G matrix --- concatenated disturbance matrix
    if sys.dist_dim > 0
        % Compute the respective extended controllability matrix (flipped)
        % [A^{N-1}F A^{N-2}F ... AF F]
        flipped_controllability_matrix_disturbance = sys.dist_mat;
        for time_index=1:time_horizon-1
            % Prepend A times (A^(time_index-1) * F) to the existing 
            % flipped_controllability_matrix_disturbance
            flipped_controllability_matrix_disturbance = ...
                [sys.state_mat * ...
                    flipped_controllability_matrix_disturbance(:, ...
                                               1:sys.dist_dim), ...
                 flipped_controllability_matrix_disturbance];
        end
        % use parts of flipped_controllability_matrix_disturbance to obtain
        %   G 
        % blocks_of_zero_required are the number of columns in the zero matrices
        %   that have state_dim number of rows. These zero matrices go to
        %   the right of the submatrices of
        %   flipped_controllability_matrix_disturbance in the G
        %   construction.
        % Start G construction from the top row
        G = [];
        for blocks_of_zero_required = time_horizon-1:-1:0
            relevant_indices_for_the_current_row = ...
                ((blocks_of_zero_required * sys.dist_dim) + 1) :... 
                (sys.dist_dim * time_horizon);
            current_row = ...
                    [ flipped_controllability_matrix_disturbance(:, ...
                          relevant_indices_for_the_current_row), ...
                      zeros(sys.state_dim, ...
                          blocks_of_zero_required * sys.dist_dim)];
            G = [G;
                 current_row];
        end
    else
        G = zeros(sys.state_dim * time_horizon, 1);
    end
end
