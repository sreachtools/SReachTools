function interp_set = interpStochReachAvoidSet(prob_thresh_of_interest,...
            polytope_above,...        
            prob_thresh_above,...
            polytope_below,...        
            prob_thresh_below)
%  Interpolate
% polytopic representations of two level sets for an (underapproximative)
% polytopic representation of a desired level set
% ============================================================================
%
% Utilize the log-concavity available at 
%
%  Abraham P. Vinod and Meeko M. K. Oishi, "TODO".
%
% Usage: TODO
%
% ============================================================================
%
% 
% Inputs:
% -------
%   prob_thresh_of_interest - Desired level set
%   polytope_above          - Polytopic representation of the higher-valued
%                               level set         
%   prob_thresh_above       - Higher valued level
%   polytope_below          - Polytopic representation of the lower-valued
%                               level set         
%   prob_thresh_below       - Lower valued level
%
% Outputs:
% --------
%   interp_set - Underapproximative polytopic representation of interpolated
%                   level set
%
% Notes:
% ------
% * NOT ACTIVELY TESTED: TODO
% * MATLAB DEPENDENCY: None
% * EXTERNAL DEPENDENCY: None
% ============================================================================
% 
% This function is part of the Stochastic Reachability Toolbox.
% License for the use of this function is given in
%      https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
%
%

    % Check they are 2 polytopes, 3 scalars in (0,1]
    inpar = inputParser();
    inpar.addRequired('prob_thresh_of_interest', @(x) validateattributes(x,...
        {'numeric'}, {'>=', 0, '<=', 1}));
    inpar.addRequired('polytope_above', @(x) validateattributes(x,...
        {'Polyhedron'}, {'nonempty'}));
    inpar.addRequired('prob_thresh_above', @(x) validateattributes(x,...
        {'numeric'}, {'>', 0, '<=', 1}));
    inpar.addRequired('polytope_thresh', @(x) validateattributes(x,...
        {'Polyhedron'}, {'nonempty'}));
    inpar.addRequired('prob_thresh_below', @(x) validateattributes(x,...
        {'numeric'}, {'>', 0, '<=', 1}));
    try
        inpar.parse(prob_thresh_of_interest,...
            polytope_above,...        
            prob_thresh_above,...
            polytope_below,...        
            prob_thresh_below);
    catch err
        exc = SrtInvalidArgsError.withFunctionName();
        exc = exc.addCause(err);
        throwAsCaller(exc);
    end

    if not(polytope_above.hasVRep && polytope_below.hasVRep)
        % Check if they have the same dimensions (needed for contains)
        exc = SrtInvalidArgsError.withFunctionName(['Polyhedra must have ',...
            'vertex representation. Please call polyhedra.computeVRep() ',...
            'before passing them as arguments.']);
        throwAsCaller(exc);
    elseif size(polytope_above.V,2) ~= size(polytope_below.V,2)
        % Check if they have the same dimensions (needed for contains)
        exc = SrtInvalidArgsError.withFunctionName(['Polyhedra must have ',...
            'same dimensions.']);
        throwAsCaller(exc);
    elseif size(polytope_above.V,1) ~= size(polytope_below.V,1)
        % Warn if polytopes do not have the same number of vertices.
        warning(['For best results, polyhedral vertices must subtend equal ',...
            'angles about their centroid.']);
    elseif not(prob_thresh_above > prob_thresh_below &&...
        polytope_below.contains(polytope_above))
        % Check if the above and the below relations are met 
        exc = SrtInvalidArgsError.withFunctionName(['Polyhedra/probability ',...
            'thresholds do not maintain the required (inclusion/lesser ',...
            'than) relation.']);
        throwAsCaller(exc);
    end

    % Get vertices
    % Get their minimum vertex representations
    polytope_below.minVRep();
    polytope_above.minVRep();
        
    below_vertices = polytope_below.V;
    above_vertices = polytope_above.V;
        
    if polytope_below.Dim == 2
        %Interpolate the vertices of the other set to get matching sets
        above_vertices_interp = vertexInterpolate(above_vertices, below_vertices);
        below_vertices_interp = vertexInterpolate(below_vertices, above_vertices);
        
        % Sort the vertices in the order of increasing angle subtended with
        % x-axis
        % Step 1: Get the centers
        below_vertices_center = mean(below_vertices);
        above_vertices_center = mean(above_vertices);
        below_vertices_interp_center = mean(below_vertices_interp);
        above_vertices_interp_center = mean(above_vertices_interp);
        % Step 2: Get the angle made by the vertex with x-axis
        below_vertices_angles = atan2(below_vertices(:,2)-below_vertices_center(2),...
                                      below_vertices(:,1)-below_vertices_center(1));
        above_vertices_angles = atan2(above_vertices(:,2)-above_vertices_center(2),...
                                      above_vertices(:,1)-above_vertices_center(1));                                  
        below_vertices_interp_angles = atan2(below_vertices_interp(:,2)-below_vertices_interp_center(2),...
                                             below_vertices_interp(:,1)-below_vertices_interp_center(1));
        above_vertices_interp_angles = atan2(above_vertices_interp(:,2)-above_vertices_interp_center(2),...
                                             above_vertices_interp(:,1)-above_vertices_interp_center(1));                                  
        % Step 3: Sort the angles
        [~,below_sort_idx] = sort(below_vertices_angles);
        [~,above_sort_idx] = sort(above_vertices_angles);
        [~,below_sort_interp_idx] = sort(below_vertices_interp_angles);
        [~,above_sort_interp_idx] = sort(above_vertices_interp_angles);
        % Step 4: Carry forward the sorting to the original listing of
        % vertices
        below_vertices = below_vertices(below_sort_idx,:);
        above_vertices = above_vertices(above_sort_idx,:);                    
        below_vertices_interp = below_vertices_interp(below_sort_interp_idx,:);
        above_vertices_interp = above_vertices_interp(above_sort_interp_idx,:);                    
    end
    
    % Compute the theta for the convex combination
    cvx_comb_theta = (log(prob_thresh_above)-log(prob_thresh_of_interest))/...
        (log(prob_thresh_above)-log(prob_thresh_below));
    
    % Use original above vertices and use interpolated below vertices
    interp_vertex_above = zeros(size(above_vertices,1), size(polytope_above.V,2));    
    for vertex_indx=1:size(above_vertices,1)
        below_vertex_interp = below_vertices_interp(vertex_indx, :);
        above_vertex = above_vertices(vertex_indx, :);
        interp_vertex_above(vertex_indx, :) = below_vertex_interp * cvx_comb_theta + ...
            (1-cvx_comb_theta) * above_vertex;
    end
    
    % Use original below vertices and use interpolated above vertices
    interp_vertex_below = zeros(size(below_vertices,1), size(polytope_above.V,2));    
    for vertex_indx=1:size(below_vertices,1)
        below_vertex = below_vertices(vertex_indx, :);
        above_vertex_interp = above_vertices_interp(vertex_indx, :);
        interp_vertex_below(vertex_indx, :) = below_vertex * cvx_comb_theta + ...
            (1-cvx_comb_theta) * above_vertex_interp;
    end
    
    % Convex hull of all vertices
    interp_set = Polyhedron('V',[interp_vertex_above;interp_vertex_below]);
end

function [v_interp] = vertexInterpolate(verts, other_verts)
    % Expects both arguments to be a collection of vertices arranged row-wise
    % verts       --- Existing vertices
    % other_verts --- Angle-generating vertices to which we interpolate

    % Angles made by the vertices of the below polytope at the origin
    center_pt = mean(verts);
    
    % Construct angles from other vertices
    angles = atan2(other_verts(:,2)-center_pt(2),other_verts(:,1)-center_pt(1));
    
    % Initalization
    v_interp = zeros(length(angles),2);
    
    % MATLAB iterator requires options to be column vectors
    fprintf('Interpolating vertices: %3d/%3d',0,length(angles));
    for angle_indx = 1:length(angles)
        % Iteration angle
        angle_iter = angles(angle_indx);
        
        %% Solve the optimization variable
        % maximize r
        % subject to center_pt + r * dir = bound_pt
        %            bound_pt in ConvexHull(vertices)
        % ConvexHull(vertices) is enforced by bound_pt = sum_i lambda_i
        % vertices_i, sum_i lambda_i = 1
        cvx_begin quiet
            variable r nonnegative;
            variable bound_pt(1,2);
            variable lambda_vec(1, length(verts)) nonnegative;
            maximize r
            subject to
                center_pt + r * [cos(angle_iter), sin(angle_iter)] == bound_pt;
                bound_pt == lambda_vec * verts;
                sum(lambda_vec) == 1;
        cvx_end        
        v_interp(angle_indx,:) = bound_pt;
        fprintf('\b\b\b\b\b\b\b%3d/%3d',angle_indx,length(angles));
    end
    fprintf('\n');
end