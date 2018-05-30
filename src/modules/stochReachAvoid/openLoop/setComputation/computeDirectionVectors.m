function set_of_direction_vectors = computeDirectionVectors(...
                                                   no_of_direction_vectors, ...
                                                   state_dimension, ...
                                                   affine_hull_of_interest_2D)
% SReachTools/stochasticReachAvoid/computeDirectionVectors: Sample a set of
% direction vectors constrained to lie within affine_hull_of_interest_2D
% =============================================================================
%
% computeDirectionVectors samples affine_hull_of_interest_2D for a set of
% "equally spaced" set of direction vectors. These vectors help in the
% construction of the polytopic underapproximation of the stochastic reach-avoid
% set in getUnderapproxStochReachAvoidSet.
%
% Any vector in the affine hull (Ae*x = be) can be described as 
%        x_0 + d_i,  i= {1,2, ...,no_of_direction_vectors}
% where,
%         Ae * v_i = 0,
%              d_i = a_i1 * v_1 + a_i2 * v_2, 
%         Ae * v_1 = Ae * v_2 = 0
% We have only v_1 and v_2 since affine_hull_of_interest_2D is assumed to be
% hyperplane that of dimension state_dimension-2.
% We can obtain v_1 and v_2 from the columns of null(Ae). 
% This function computes d_i. 
%
% USAGE: See stochasticReachAvoid/getUnderapproxStochReachAvoidSet
%
% =============================================================================
%
% set_of_direction_vectors = computeDirectionVectors(...
%                                                 no_of_direction_vectors, ...
%                                                 state_dimension, ...
%                                                 affine_hull_of_interest_2D)
% 
% Inputs:
% -------
%   no_of_direction_vectors    - Number of unique directions defining the polytope
%                                vertices 
%   state_dimension            - Dimension of the state space
%                                (For example, sys.state_dimension)
%   affine_hull_of_interest_2D - Affine hull whose slice of the stochastic
%                                reach-avoid set is of interest
%                                Dimension state_dimension-2
%                                Define this by Polyhedron('He',[A_eq, b_eq])
%
% Outputs:
% --------
%   set_of_direction_vectors   - Set of direction vectors sampled from
%                                affine_hull_of_interest_2D
%
% See also getUnderapproxStochReachAvoidSet
%
% Notes:
% ------
% * NOT ACTIVELY TESTED: Builds on other tested functions.
% * EXTERNAL DEPENDENCY: Uses MPT3
%                        Needs MPT3 for affine_hull_of_interest_2D manipulation
% =============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/abyvinod/SReachTools/blob/master/LICENSE
%
%

    % Check no_of_direction_vectors is a positive scalar
    assert( isscalar(no_of_direction_vectors) &&...
                no_of_direction_vectors > 0, ...
            'SReachTools:invalidArgs', ...
            'No_of_direction_vectors needs to be a positive scalar');

    % Check state_dimension is a positive scalar
    assert( isscalar(state_dimension) &&...
                state_dimension > 0, ...
            'SReachTools:invalidArgs', ...
            'State dimension needs to be a positive scalar');

    % Check affine_hull_of_interest_2D_A has rank 2, has sys.state_dimension
    % columns, and is non-empty
    affine_hull_of_interest_2D.minAffineRep();
    assert( rank(affine_hull_of_interest_2D.Ae) == state_dimension - 2, ...
            'SReachTools:invalidArgs', ...
            'Dimension of affine hull is not sys.state_dimension-2');
    assert( affine_hull_of_interest_2D.Dim == state_dimension, ...
            'SReachTools:invalidArgs', ...
            'Affine hull is not embedded in n-dimensional space');
    assert(~affine_hull_of_interest_2D.isEmptySet(), ...
            'SReachTools:invalidArgs', ...
            'Affine hull is empty');

    % Compute the direction angles (theta_vector)
    theta_vector_extra_elem = linspace(0, 2*pi, no_of_direction_vectors + 1);
    theta_vector = theta_vector_extra_elem(1:end-1);
    % Compute the direction vectors (cos(theta), sin(theta)) corresponding to
    % theta_vector
    coefficients_for_basis = [cos(theta_vector); 
                              sin(theta_vector)];
    % Compute two vectors spanning null space of affine_hull_of_interest_2D.Ae
    basis_vectors_for_null_space = null(affine_hull_of_interest_2D.Ae);

    % d_i = [v_1 v_2]*[cos(theta_1) cos(theta_2) cos(theta_{no_of_d_i});
    %                  sin(theta_1) sin(theta_2) sin(theta_{no_of_d_i})];
    set_of_direction_vectors = basis_vectors_for_null_space * ...
                                                coefficients_for_basis;
    %TODO: Need to check for affine_hull_of_interest_2D w/ non-Cartesian planes
end

% LINPROG-based check for feasibility
% linprogOptions = optimoptions('linprog','Display','off');
% [~,~,exit_flag_for_linprog]=linprog(zeros(state_dimension,1), ...
%                                     [],[], ...
%                                     affine_hull_of_interest_2D.Ae, ...
%                                     affine_hull_of_interest_2D.be, ...
%                                     [],[], ...
%                                     linprogOptions);
