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
        warning('For best results, polyhedra must have same vertices.');
        %% TODO: Do they need to extrapolate to xmax?
    elseif not(prob_thresh_above > prob_thresh_below &&...
        polytope_below.contains(polytope_above))
        keyboard
        % Check if the above and the below relations are met 
        exc = SrtInvalidArgsError.withFunctionName(['Polyhedra/probability ',...
            'thresholds do not maintain the required (inclusion/lesser ',...
            'than) relation.']);
        throwAsCaller(exc);
    end

    %% Compute the theta for the convex combination
    cvx_comb_theta = (log(prob_thresh_above)-log(prob_thresh_of_interest))/...
        (log(prob_thresh_above)-log(prob_thresh_below));

    %% Interpolate each vertex point based on the above theta rule
    n_interp_vertex = min(size(polytope_above.V,1),size(polytope_below.V,1));
    interp_vertex = zeros(n_interp_vertex, size(polytope_above.V,2));
    for vertex_indx=1:n_interp_vertex
        below_vertex = polytope_below.V(vertex_indx, :);
        above_vertex = polytope_above.V(vertex_indx, :);
        interp_vertex(vertex_indx, :) = below_vertex * cvx_comb_theta + ...
            (1-cvx_comb_theta) * above_vertex;
    end
    interp_set = Polyhedron('V',interp_vertex);
end
