%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Requires: Ellipsoidal toolbox
%
% Tests: 1) get_FSRPD_mean_covariance_matrix, and
%        2) get_ctrb_and_state_transition_matrices_unicycle, and
%
% Description: Computes the mean trajectory and 99% confidence ellipses for 
%       a unicycle whose velocity is taken from a Gaussian random vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
disp('Note that ellipses here are just a representation of \Sigma and are not avoid sets!')
disp('Produces Figure 1 (Infinity symbol with increasing ellipse volumes)');
% Disable verbosity of ellipsoidal toolbox
evalc('ellipsoids_init');
global ellOptions;
ellOptions.verbose = 0;
% Compute the ctrb and state transition matrices based on unicycle dynamics
mode_sequence_omega = [0,0,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,0,0,0,0,0,0,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5];
conf_level = 0.99;
% mode_sequence_omega = zeros(1,5);
theta_init = pi/4;
sampling_time = 0.05;
% Random velocity parameters
mu_velocity = 5;
variance_velocity = 1;
sys = getDubinsCarLtv('vel-dist',...
    mode_sequence_omega,...
    theta_init,...
    sampling_time,...
    RandomVector('Gaussian', mu_velocity, variance_velocity));

% Initial location
x_init = 0;
y_init = 0;
obstacle_init_location = [x_init;y_init];
% obstacle_init_location = RandomVector('Gaussian', zeros(2,1), 1e-2*eye(2));

% Plot axes
xmin=-3.5;
xmax=3;
ymin=-2.5;
ymax=1.5;

elapsed_time = zeros(1, length(mode_sequence_omega)+1);

figure(1)
clf
hold on
% Note that time flows from 0 to time_horizon=length(mode_sequence_omega)
% or in MATLAB parlance, 1 to length(mode_sequence_omega)+1
for time_of_interest = 1:length(mode_sequence_omega)+1
    
    timerVal=tic;
    elapsed_time(time_of_interest)=toc(timerVal);
    [mu_trajectory{time_of_interest}, sigma_obstacle{time_of_interest}] =...
        getFSRPDMeanCov(sys, obstacle_init_location, time_of_interest);
    
    % Plot the mean trajectory
    scatter(mu_trajectory{time_of_interest}(1),mu_trajectory{time_of_interest}(2),50,'bx');
    if ~issymmetric(sigma_obstacle{time_of_interest})
        % Sometimes numerical imprecision to the order of 10^(-19) happens
        if abs(sigma_obstacle{time_of_interest}(2)-sigma_obstacle{time_of_interest}(3)) < eps
            % Neutralize it
            sigma_obstacle{time_of_interest} = (sigma_obstacle{time_of_interest} + sigma_obstacle{time_of_interest}')/2;        
        else
            % Should never happen! Something is wrong
            error('The sigma_obstacle matrix is not symmetric');
        end
    end
    % 99% confidence ellipses creation. See CDC2017 Gleason et. al, Section
    % IV.A. Equations (32) and (33)
    try
        E = ellipsoid(mu_trajectory{time_of_interest},sigma_obstacle{time_of_interest}*chi2inv(conf_level,2));
        plot(E);
    catch
        E = ellipsoid(mu_trajectory{time_of_interest},zeros(2,2));
        plot(E);
    end
end
set(gca,'FontSize',20);
%title(sprintf('Mean and Covariance matrices\n(99%% confidence ellipses)'));
fsprintf('Mean and Covariance matrices\n(99%% confidence ellipses)\n');
% Blue circles mark the mean\nRed ellipses mark 99%% confidence ellipses'));
leg=legend('Mean','99% confidence ellipse');
grid on
axis equal
axis([xmin xmax ymin ymax])
set(gca,'XTick',[xmin:0.5:xmax])
set(gca,'YTick',[ymin:0.5:ymax])
set(leg,'Location','SouthEast')
% Check if circle
%plot(ellipsoid([0;0],0.00100*eye(2)));
box on
%figure(106)
%clf
%stem(elapsed_time)
