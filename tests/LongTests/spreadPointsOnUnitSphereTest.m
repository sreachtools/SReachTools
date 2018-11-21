%% Check arrangment via difference of convex programming    
% clear;close all;clc;cvx_clear;

%% Check for a particular arrangement
verbose = 1;
n_dim = 3;
n_points = 2^(n_dim) * 10 + 2*n_dim;

%tic
[opt_locations, separation]=spreadPointsOnUnitSphere(n_dim, n_points, verbose);
%toc

%% Uncomment these lines if you wish to plot
if n_dim == 2 || n_dim == 3
    figure();
    if n_dim == 2
        quiver(zeros(1,n_points), zeros(1,n_points),...
            opt_locations(1,:),opt_locations(2,:));
    else
        quiver3(zeros(1,n_points),zeros(1,n_points),zeros(1,n_points),...
            opt_locations(1,:),opt_locations(2,:),opt_locations(3,:));
    end
    axis equal;
    grid on;
    box on;
end
