function poly_array = getDynProgLevelSets2D(cell_of_grid_x, prob_x, ...
    prob_lvls, safety_tube)
% SReachTools/stochasticReachAvoid/getDynProgLevelSets2D Get level sets based 
% on the value function returned by getDynProgSolForTube
% ============================================================================
%
% The function computes an array of polytopes based on the results from
% getDynProgSolForTube
%
% See also examples/doubleIntegratorDynamicProgramming.m, SReachDynProg
%
% ============================================================================
%
% poly_array = getDynProgLevelSets2D(sys, prob_x, prob_lvls)
% 
% Inputs:
% -------
%   cell_xvec  - Gridding along the particular dimension (sys.state_dim x 1
%                cell array, with grid info along each dimension) | Output for
%                SReachDynProg
%   prob_x     - Probability values at each grid point (M number of them) in
%                grid_X (Mx1 array)
%   prob_lvls  - A vector containing safety probability thresholds of interest
%                Each element needs to be within [0,1].
%   safety_tube- Safety tube used for the dynamic programming solution
%
% Outputs:
% --------
%   poly_array - Array of Polyhedron objects for each of the level sets
%
% Notes:
% ------
% * Uses MATLAB's contour function to create the boundary info which is then fed
%   to MPT (takes a convex hull). Because contour function may ignore corners
%   (if value function saturates), we check for all the corners when
%   constructing the polytope.
% * To be used in conjunction with SReachDynProg
% 
% ============================================================================
% 
%   This function is part of the Stochastic Reachability Toolbox.
%   License for the use of this function is given in
%        https://sreachtools.github.io/license/
%
%

    %% Input handling
    validateattributes(cell_of_grid_x, {'cell'},{'vector'})
    validateattributes(prob_x, {'numeric'}, {'vector','nonempty'})
    validateattributes(prob_lvls, {'numeric'}, {'vector','nonempty'})
    validateattributes(safety_tube, {'Tube'}, {'nonempty'});
    % Check if prob_lvls is a [0,1] vector (check min >= 0 and max <=1)
    if min(prob_lvls) <0 || max(prob_lvls) > 1
        throwAsCaller(SrtInvalidArgsError('prob_lvls need to be [0,1].'));    
    end
    % Check if cell_of_grid_x has two elements
    if length(cell_of_grid_x)~=2
        throwAsCaller(SrtInvalidArgsError(['getDynProgLevelSets is meant', ...
            'only for 2-D systems (3D cell_of_grid_x given)!']));
    end

    %% Computation
    % x grid points
    x1vec = cell_of_grid_x{1};
    x2vec = cell_of_grid_x{2};
    
    % Corners of target set at t=0
    corners = safety_tube(1).outerApprox.V;
    
    % Initialize the set of vertices for the polytopes
    poly_array_vertices = cell(length(prob_lvls),1);

    % Obtain the matrix form of prob_x
    grid_probability_mat = reshape(prob_x, length(x2vec),[]);

    % Obtain the contour plot
    if length(prob_lvls)==1
        % Matlab's contourc requires [p p] if only one level is of interest
        [C_DP]=contourc(x1vec, x2vec, grid_probability_mat, ...
            [prob_lvls prob_lvls]);
    else
        [C_DP]=contourc(x1vec, x2vec, grid_probability_mat, prob_lvls);
    end

    % Parse the contour matrix
    col_indx = 1;
    while col_indx <= length(C_DP)
        % Contour matrix has a specific structure. See
        % https://www.mathworks.com/help/matlab/ref/matlab.graphics.chart.primitive.contour-properties.html#d119e148367_panel-group
        current_level = C_DP(1,col_indx);
        current_level_indx = find(abs(prob_lvls - current_level)<eps);
        no_of_points = C_DP(2,col_indx);
        poly_array_vertices{current_level_indx} =...
            [poly_array_vertices{current_level_indx}, ...
             C_DP(:,col_indx+1:col_indx+no_of_points)];
        col_indx = col_indx + no_of_points + 1;        
    end

    % Construct the polyhedrons
    poly_array = repmat(Polyhedron.emptySet(2),length(prob_lvls),1);
    for poly_indx = 1:length(prob_lvls)
        originalPolytope = Polyhedron('V', poly_array_vertices{poly_indx}');
        originalPolytope.minVRep();
        originalVertices = originalPolytope.V;
        n_orig_vertices  = size(originalVertices,1);            
        for corner_indx = 1:size(corners,1)
            tempPolytope = Polyhedron('V',[originalVertices;
                                           corners(corner_indx,:)]);
            tempPolytope.minVRep();
            if (size(tempPolytope.V,1) == n_orig_vertices + 1) ||...
                (size(tempPolytope.V,1) == n_orig_vertices)
                warning('SReachTools:runtime', sprintf(['MATLAB''s ', ...
                    'contour matrix missed a corner! \nAdding (%d, %d) to ', ...
                    'polytope vertex list for level=%1.2f.'], ...
                    corners(corner_indx,:)', prob_lvls(poly_indx)));
                poly_array_vertices{poly_indx} = ...
                    [poly_array_vertices{poly_indx}, corners(corner_indx,:)'];
            end
        end
        poly_array(poly_indx) = Polyhedron('V',poly_array_vertices{poly_indx}');        
    end
end
