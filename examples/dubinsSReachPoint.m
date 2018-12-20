%% Controller synthesis using |SReachPoint| for a Dubin's vehicle
% This example will demonstrate the use of |SReachTools| for controller
% synthesis in a stochastic continuous-state discrete-time linear time-varying
% (LTV) systems. This example script is part of the |SReachTools| toolbox, which
% is licensed under GPL v3 or (at your option) any later version. A copy of this
% license is given in
% <https://github.com/unm-hscl/SReachTools/blob/master/LICENSE
% https://github.com/unm-hscl/SReachTools/blob/master/LICENSE>.
% 
% In this example script, we discuss how to use |SReachPoint| to synthesize
% open-loop controllers and affine-disturbance feedback controllers for the
% problem of stochastic reachability of a target tube. We demonstrate the
% following solution techniques:
% 
% * |chance-open|: Chance-constrained approach that uses risk allocation and 
%    piecewise-affine approximations to formulate a linear program to
%    synthesize an open-loop controller (See
%    <http://hscl.unm.edu/affinecontrollersynthesis Vinod and Oishi, Hybrid
%    Systems: Computation and Control, 2019 (submitted)>,
%    <http://doi.org/10.1109/CDC.2013.6760626 Lesser et. al., Conference on
%    Decision and Control, 2013>)
% * |genzps-open|: Fourier transforms that uses
%    <http://www.math.wsu.edu/faculty/genz/software/matlab/qsimvnv.m Genz's 
%    algorithm> to formulate a nonlinear log-concave optimization problem to be
%    solved using MATLAB's patternsearch to synthesize an open-loop controller
%    (See <http://doi.org/10.1109/LCSYS.2017.2716364 Vinod and Oishi, Control
%    System Society- Letters, 2017>)
% * |particle-open|: Particle control filter approach that formulates a 
%    mixed-integer linear program to synthesize an open-loop controller (See
%    <http://doi.org/10.1109/CDC.2013.6760626 Lesser et. al., Conference on
%    Decision and Control, 2013>)
% * |voronoi-open|: Particle control filter approach that formulates a
%    mixed-integer linear program to synthesize an open-loop controller. In
%    contrast to |particle-open|, |voronoi-open| permits a user-specified upper
%    bound on the overapproximation error in the maximal reach probability and
%    has significant computational advantages due to its undersampling approach.
%    (See <arxiv_link_TODO Sartipizadeh et. al., American Control Conference,
%    2019 (submitted)>)
% * |chance-affine|: Chance-constrained approach that uses risk allocation and 
%    piecewise-affine approximations to formulate a difference-of-convex program
%    to synthesize a closed-loop (affine disturbance feedback) controller.  The
%    controller synthesis is done by solving a series of second-order cone
%    programs. (See <http://hscl.unm.edu/affinecontrollersynthesis Vinod and
%    Oishi, Hybrid Systems: Computation and Control, 2019 (submitted)>)
%
% All computations were performed using MATLAB on an Intel Xeon CPU with 3.4GHz
% clock rate and 32 GB RAM. The simulation times for individual methods are
% reported in each section along with a Monte-Carlo simulation validation. The
% overall simulation time was 16 minutes. For sake of clarity, all commands were 
% asked to be verbose (via SReachPointOptions). In practice, this can be turned 
% off.
%

% Prescript running: Initializing srtinit, if it already hasn't been initialized
close all;clearvars;srtinit;srtinit --version;

%% Problem formulation: Stochastic reachability of a target tube
% Given an initial state $x_0$, a time horizon $N$, a linear system dynamics
% $x_{k+1} = A_k x_k + B_k u_k + F w_k$ for $k\in \{0,1,...,N-1\}$, and a target
% tube ${\{\mathcal{T}_k\}}_{k=0}^N$, we wish to design an admissible controller
% that maximizes the probability of the state staying with the target tube. This
% maximal reach probability, denoted by $V^\ast(x_0)$, is obtained by solving
% the following optimization problem
%
% $$ V^\ast(x_0) = \max_{\overline{U}\in \mathcal{U}^N} P^{x_0,
% \overline{U}}_{X} \{ \forall k, x_k \in \mathcal{T}_k\}.$$
%
% Here, $\overline{U}$ refers to the control policy which satisfies the control
% bounds specified by the input space $\mathcal{U}$ over the entire time
% horizon $N$, $X= {[x_1\ x_2\ \ldots\ x_N]}^\top$ is the concatenated state
% vector, and the target tube is a sequence of sets
% ${\{\mathcal{T}_k\}}_{k=0}^N$.  Here, $X$ is a random vector with probability
% measure $P^{x_0,\overline{U}}_X$ which is a parameterized by the initial state
% $x_0$ and policy $\overline{U}$.  
%
% In the general formulation requires $\overline{U}$ is given by a sequence
% of (potentially time-varying and nonlinear) state-feedback controllers. To
% compute such a policy, we have to resort to dynamic programming which suffers
% from the curse of dimensionality. See these papers for details
% <https://doi.org/10.1016/j.automatica.2008.03.027 Abate et. al, Automatica,
% 2008>, <https://doi.org/10.1016/j.automatica.2010.08.006 Summers and Lygeros,
% Automatica, 2010>, and <https://arxiv.org/abs/1810.05217 Vinod and Oishi,
% IEEE Trans. Automatic Control, 2018 (submitted)>.
%
% |SReachPoint| provides multiple ways to compute an *underapproximation* of
% $V^\ast(x_0)$ by restricting the search to the following controllers:
% 
% * open-loop controller: The controller provides a sequence of control actions
% $\overline{U}={[u_0\ u_1\ \ldots\ u_{N-1}]}^\top\in \mathcal{U}^N$
% parameterized only by the initial state. This controller does not account for
% the actual state realization and therefore can be conservative. However,
% computing this control sequence is easy due to known convexity properties of
% the problem. See <https://arxiv.org/abs/1810.05217 Vinod and Oishi,
% IEEE Trans. Automatic Control, 2018 (submitted)> for more details.  Apart from
% |particle-open|, all approaches provide guaranteed underapproximations or
% underapproximations to a user-specifed error.
% * affine-disturbance feedback controller: The controller is a characterized by
% an affine transformation of the concatenated disturbance vector. The gain
% matrix is forced to be lower-triangular for the causality, resulting in the
% control action at $k$ be dependent only the past disturbance values. 
% Here, the control action at time $k\in \{0,1,\ldots,N-1\}$ is given by $u_k =
% \sum_{i=0}^{k-1} M_{ki} w_i + d_k$.  We optimize for $M_{ki}$ and $d_k$ for
% every $k,i$, and the controller is given by $\overline{U}=MW + d\in
% \mathcal{U}^N$, with $W=[w_0\ w_1\ \ldots\ w_{N-1}]$ denoting the concatenated
% disturbance random vector. By construction, $\overline{U}$ is now random, and
% it can not satisfy hard control bounds with non-zero $M_{ki}$ and unbounded
% $W$. Therefore, we relax the control bound constraints
% $\overline{U}\in\mathcal{U}^N$ to a chance constraint, $P_W\{MW + d\in
% \mathcal{U}^N\}\geq 1-\Delta_U$ permitting the user to specify the
% probabilistic violation $\Delta_U\in[0,1)$ of the control bounds. We then
% construct a lower bound for the maximal reach probability when the affine
% disturbance feedback controller is used under saturation to meet the hard
% control bounds. In contrast to the open-loop controller synthesis, affine
% disturbance feedback controller synthesis is a non-convex problem, and we
% obtain a locally optimal solution using difference-of-convex programming. 
% See <http://hscl.unm.edu/affinecontrollersynthesis
% Vinod and Oishi, Hybrid Systems: Computation and Control, 2019 (submitted)>
% for more details.
%
% All of our approaches are grid-free resulting in highly scalable solutions,
% especially for Gaussian-perturbed linear systems. 
%
% In this example, we perform controller synthesis that maximizes the
% probability of a Dubin's vehicle to stay within a time-varying collection of
% target sets. We model the Dubin's vehicle with known turning rate sequence as
% a linear time-varying system.

%% Dubin's vehicle dynamics
% We consider a Dubin's vehicle with known turning rate sequence
% $\overline{\omega} = {[\omega_0\ \omega_1\ \ldots\ \omega_{T-1}]}^\top
% \in R^T$, with additive Gaussian disturbance. The resulting dynamics are,
% 
% $$x_{k+1} = x_k + T_s \cos\left(\theta_0 + \sum_{i=1}^{k-1} 
% \omega_i T_s\right) v_k + \eta^x_k$$
%
% $$y_{k+1} = y_k + T_s \sin\left(\theta_0 + \sum_{i=1}^{k-1} 
% \omega_i T_s\right) v_k + \eta^y_k$$
%
% where $x,y$ are the positions (state) of the Dubin's vehicle in $\mathrm{x}$- 
% and $\mathrm{y}$- axes, $v_k$ is the velocity of the vehicle (input), 
% $\eta^{(\cdot)}_k$ is the additive Gaussian disturbance affecting the
% dynamics, $T_s$ is the sampling time, and $\theta_0$ is the initial heading
% direction.  We define the disturbance as ${[\eta^x_k\ \eta^y_k]}^\top\sim
% \mathcal{N}({[0\ 0]}^\top, 10^{-3}I_2)$.

n_mcarlo_sims = 1e5;                        % Monte-Carlo simulation particles
n_mcarlo_sims_affine = 1e5;                 % For affine controllers
sampling_time = 0.1;                        % Sampling time
init_heading = pi/10;                       % Initial heading 
% Known turning rate sequence
time_horizon = 12;
omega = pi/time_horizon/sampling_time;
half_time_horizon = round(time_horizon/2);
turning_rate = [omega*ones(half_time_horizon,1);
               -omega*ones(half_time_horizon,1)];
% Input space definition
umax = 10;
input_space = Polyhedron('lb',0,'ub',umax);
% Disturbance matrix and random vector definition
dist_matrix = eye(2);
eta_dist_gauss = RandomVector('Gaussian',zeros(2,1), 0.001 * eye(2));

sys_gauss = getDubinsCarLtv('add-dist', turning_rate, init_heading, ...
    sampling_time, input_space, dist_matrix, eta_dist_gauss);


%% Target tube definition
% We define the target tube to be a collection of time-varying boxes
% $\{\mathcal{T}_k\}_{k=0}^N$ where $N$ is the time horizon.
%
% In this problem, we define $\mathcal{T}_k$ to be centered about the nominal
% trajectory with fixed velocity of $u_\mathrm{max} * 3/2$ (faster than the
% maximum velocity allowed) and the heading angle sequence with $\pi/2$ removed.
% The half-length of these boxes decay exponentially with a time constant which
% is $N/2$.

v_nominal = umax * 2/3;                 % Nominal trajectory's heading velocity
                                        % TODO: Gaussian was 3/2
% Construct the nominal trajectory
[~,H,~] = sys_gauss.getConcatMats(time_horizon);
center_box_X = [zeros(2,1);H * (v_nominal * ones(time_horizon,1))];
center_box = reshape(center_box_X,2,[]);
% Box sizes
box_halflength_at_0 = 4;                % Box half-length at t=0
time_const = 1/2*time_horizon;          % Time constant characterize the
                                        % exponentially decaying box half-length

% Target tube definition as well as plotting
target_tube_cell = cell(time_horizon + 1,1); % Vector to store target sets
figure(100);clf;hold on
for itt = 0:time_horizon
    % Define the target set at time itt
    target_tube_cell{itt+1} = Polyhedron(...
        'lb',center_box(:, itt+1) -box_halflength_at_0*exp(- itt/time_const),...
        'ub', center_box(:, itt+1) + box_halflength_at_0*exp(- itt/time_const));
    if itt==0
        % Remember the first the tube
        h_target_tube = plot(target_tube_cell{1},'alpha',0.5,'color','y');
    else
        plot(target_tube_cell{itt+1},'alpha',0.08,'LineStyle',':','color','y');
    end            
end
axis equal        
h_nominal_traj = scatter(center_box(1,:), center_box(2,:), 50,'ks','filled');        
h_vec = [h_target_tube, h_nominal_traj];
legend_cell = {'Target tube', 'Nominal trajectory'};
legend(h_vec, legend_cell, 'Location','EastOutside', 'interpreter','latex');
xlabel('x');
ylabel('y');
axis equal
box on;
grid on;
drawnow;

% Target tube definition
target_tube = Tube(target_tube_cell{:});

%% Specifying initial states and which options to run
chance_open_run_gauss = 0;
genzps_open_run_gauss = 0;
particle_open_run_gauss = 0;
voronoi_open_run_gauss = 0;
chance_affine_run_gauss = 0;

% Initial states for each of the method
init_state_chance_open_gauss = [2;2] + [-1;-1];
init_state_genzps_open_gauss = [2;2] + [1;-1];
init_state_particle_open_gauss = [2;2] + [0;1];
init_state_voronoi_open_gauss = [2;2] + [1.5;1.5];
init_state_chance_affine_gauss = [2;2] + [2;1];


%% Quantities needed to compute the optimal mean trajectory 
% We first compute the dynamics of the concatenated state vector $X = Z x_0
% + H U + G W$, and compute the concatentated random vector $W$ and its mean.
[Z,H,G] = sys_gauss.getConcatMats(time_horizon);
% Compute the mean trajectory of the concatenated disturbance vector
muW_gauss = sys_gauss.dist.concat(time_horizon).mean();

%% |SReachPoint|: |chance-open|
% This method is discussed in <http://hscl.unm.edu/affinecontrollersynthesis
% Vinod and Oishi, Hybrid Systems: Computation and Control, 2019 (submitted)>.
% It was introduced for stochastic reachability in
% <http://doi.org/10.1109/CDC.2013.6760626 Lesser et. al., Conference on
% Decision and Control, 2013>.
%
% This approach implements the chance-constrained approach to compute an optimal 
% open-loop controller. It uses risk allocation and piecewise-affine
% overapproximation of the inverse normal cumulative density function to
% formulate a linear program for this purpose. Naturally, this is one of the
% fastest ways to compute an open-loop controller and an underapproximative
% probabilistic guarantee of safety. However, due to the use of Boole's
% inequality for risk allocation, it provides a conservative estimate of safety
% using the open-loop controller.

if chance_open_run_gauss
    fprintf('\n\nSReachPoint with chance-open\n');
    % Set the maximum piecewise-affine overapproximation error to 1e-3
    opts = SReachPointOptions('term', 'chance-open','pwa_accuracy',1e-3);
    timerVal=tic;
    [prob_chance_open_gauss, opt_input_vec_chance_open_gauss] = SReachPoint('term', ...
        'chance-open', sys_gauss, init_state_chance_open_gauss, target_tube, opts);
    elapsed_time_chance_open_gauss = toc(timerVal);
    if prob_chance_open_gauss
        % Optimal mean trajectory construction
        % mean_X = Z * x_0 + H * U + G * \mu_W
        opt_mean_X_chance_open_gauss = Z * init_state_chance_open_gauss + ...
            H * opt_input_vec_chance_open_gauss + G * muW_gauss;
        opt_mean_traj_chance_open_gauss = reshape(opt_mean_X_chance_open_gauss, ...
            sys_gauss.state_dim,[]);
        % Check via Monte-Carlo simulation
        concat_state_realization = generateMonteCarloSims(n_mcarlo_sims, ...
            sys_gauss, init_state_chance_open_gauss, time_horizon,...
            opt_input_vec_chance_open_gauss);
        mcarlo_result = target_tube.contains(concat_state_realization);
        simulated_prob_chance_open_gauss = sum(mcarlo_result)/n_mcarlo_sims;
    else
        simulated_prob_chance_open_gauss = NaN;
    end
    fprintf('SReachPoint underapprox. prob: %1.2f | Simulated prob: %1.2f\n',...
        prob_chance_open_gauss, simulated_prob_chance_open_gauss);
    fprintf('Computation time: %1.3f\n', elapsed_time_chance_open_gauss);
end

%% |SReachPoint|: |genzps-open|
% This method is discussed in <http://doi.org/10.1109/LCSYS.2017.2716364
% Vinod and Oishi, Control System Society- Letters, 2017>.
%
% This approach implements the Fourier transform-based approach to compute an
% optimal open-loop controller. It uses
% <http://www.math.wsu.edu/faculty/genz/software/matlab/qsimvnv.m Genz's
% algorithm> to compute the probability of safety and optimizes the joint chance
% constraint involved in maximizing this probability. To handle the noisy
% behaviour of the Genz's algorithm, we rely on MATLAB's |patternsearch| for the
% nonlinear optimization. Internally, we use the
% |chance-open| to initialize the nonlinear solver. Hence, this approach will
% return an open-loop controller with safety at least as good as |chance-open|.
if genzps_open_run_gauss
    fprintf('\n\nSReachPoint with genzps-open\n');
    opts = SReachPointOptions('term', 'genzps-open', ...
        'PSoptions',psoptimset('display','iter'));
    timerVal = tic;
    [prob_genzps_open_gauss, opt_input_vec_genzps_open_gauss] = SReachPoint('term', ...
        'genzps-open', sys_gauss, init_state_genzps_open_gauss, target_tube, opts);
    elapsed_time_genzps_gauss = toc(timerVal);
    if prob_genzps_open_gauss > 0
        % Optimal mean trajectory construction
        % mean_X = Z * x_0 + H * U + G * \mu_W
        opt_mean_X_genzps_open_gauss =  Z * init_state_genzps_open_gauss + ...
            H * opt_input_vec_genzps_open_gauss + G * muW_gauss;
        opt_mean_traj_genzps_open_gauss = reshape(opt_mean_X_genzps_open_gauss, ...
            sys_gauss.state_dim,[]);
        % Check via Monte-Carlo simulation
        concat_state_realization = generateMonteCarloSims(n_mcarlo_sims, ...
            sys_gauss, init_state_genzps_open_gauss, time_horizon,...
            opt_input_vec_genzps_open_gauss);
        mcarlo_result = target_tube.contains(concat_state_realization);
        simulated_prob_genzps_open_gauss = sum(mcarlo_result)/n_mcarlo_sims;
    else
        simulated_prob_genzps_open_gauss = NaN;
    end
    fprintf('SReachPoint underapprox. prob: %1.2f | Simulated prob: %1.2f\n',...
        prob_genzps_open_gauss, simulated_prob_genzps_open_gauss);
    fprintf('Computation time: %1.3f\n', elapsed_time_genzps_gauss);    
end

%% |SReachPoint|: |particle-open|
% This method is discussed in <http://doi.org/10.1109/CDC.2013.6760626
% Lesser et. al., Conference on Decision and Control, 2013>.
%
% This approach implements the particle control approach to compute an open-loop
% controller. It is a sampling-based technique and hence the resulting
% probability estimate is random with its variance going to zero as the number
% of samples considered goes to infinity. Note that since a mixed-integer linear
% program is solved underneath with the number of binary variables corresponding
% to the number of particles, using too many particles can cause an exponential
% increase in computational time.
if particle_open_run_gauss
    fprintf('\n\nSReachPoint with particle-open\n');
    opts = SReachPointOptions('term','particle-open','verbose',1,...
        'n_particles',50);
    timerVal = tic;
    [prob_particle_open_gauss, opt_input_vec_particle_open_gauss] = SReachPoint('term', ...
        'particle-open', sys_gauss, init_state_particle_open_gauss, target_tube, opts);
    elapsed_time_particle_gauss = toc(timerVal);
    if prob_particle_open_gauss > 0
        % Optimal mean trajectory construction
        % mean_X = Z * x_0 + H * U + G * \mu_W
        opt_mean_X_particle_open_gauss =  Z * init_state_particle_open_gauss + ...
            H * opt_input_vec_particle_open_gauss + G * muW_gauss;
        opt_mean_traj_particle_open_gauss =...
            reshape(opt_mean_X_particle_open_gauss, sys_gauss.state_dim,[]);
        % Check via Monte-Carlo simulation
        concat_state_realization = generateMonteCarloSims(n_mcarlo_sims, ...
            sys_gauss, init_state_particle_open_gauss, time_horizon,...
            opt_input_vec_particle_open_gauss);
        mcarlo_result = target_tube.contains(concat_state_realization);
        simulated_prob_particle_open_gauss = sum(mcarlo_result)/n_mcarlo_sims;
    else
        simulated_prob_particle_open_gauss = NaN;
    end
    fprintf('SReachPoint approx. prob: %1.2f | Simulated prob: %1.2f\n',...
        prob_particle_open_gauss, simulated_prob_particle_open_gauss);
    fprintf('Computation time: %1.3f\n', elapsed_time_particle_gauss);
end

%% |SReachPoint|: |voronoi-open|
% This method is discussed in <arxiv_link_TODO Sartipizadeh et. al., 
% American Control Conference, 2019 (submitted)>
%
% This approach implements the undersampled particle control approach to compute
% an open-loop controller. It computes, using k-means, a representative sample
% realization of the disturbance which is significantly smaller. This
% drastically improves the computational efficiency of the particle control
% approach. Further, because it uses Hoeffding's inequality, the user can
% specify an upper-bound on the overapproximation error. The undersampled
% probability estimate is used to create a lower bound of the solution
% corresponding to the original particle control problem with appropriate
% (typically large) number of particles. Thus, this has all the benefits of the
% |particle-open| option, with additional benefits of being able to specify a
% maximum overapproximation error as well being computationally tractable.
if voronoi_open_run_gauss
    fprintf('\n\nSReachPoint with voronoi-open\n');
    opts = SReachPointOptions('term','voronoi-open','verbose',1,...
        'max_overapprox_err', 1e-2);
    timerVal = tic;
    [prob_voronoi_open_gauss, opt_input_vec_voronoi_open_gauss] = SReachPoint('term', ...
        'voronoi-open', sys_gauss, init_state_voronoi_open_gauss, target_tube, opts);
    elapsed_time_voronoi_gauss = toc(timerVal);
    if prob_voronoi_open_gauss > 0
        % Optimal mean trajectory construction
        % mean_X = Z * x_0 + H * U + G * \mu_W
        opt_mean_X_voronoi_open_gauss =  Z * init_state_voronoi_open_gauss + ...
            H * opt_input_vec_voronoi_open_gauss + G * muW_gauss;
        opt_mean_traj_voronoi_open_gauss =...
            reshape(opt_mean_X_voronoi_open_gauss, sys_gauss.state_dim,[]);
        % Check via Monte-Carlo simulation
        concat_state_realization = generateMonteCarloSims(n_mcarlo_sims, ...
            sys_gauss, init_state_voronoi_open_gauss,time_horizon,...
            opt_input_vec_voronoi_open_gauss);
        mcarlo_result = target_tube.contains(concat_state_realization);
        simulated_prob_voronoi_open_gauss = sum(mcarlo_result)/n_mcarlo_sims;
    else
        simulated_prob_voronoi_open_gauss = NaN;
    end
    fprintf('SReachPoint approx. prob: %1.2f | Simulated prob: %1.2f\n',...
        prob_voronoi_open_gauss, simulated_prob_voronoi_open_gauss);
    fprintf('Computation time: %1.3f\n', elapsed_time_voronoi_gauss);
end

%% |SReachPoint|: |chance-affine|
% This method is discussed in <http://hscl.unm.edu/affinecontrollersynthesis
% Vinod and Oishi, Hybrid Systems: Computation and Control, 2019 (submitted)>.
%
% This approach implements the chance-constrained approach to compute a locally
% optimal affine disturbance feedback controller. In contrast to |chance-open|,
% this approach optimizes for an affine feedback gain for the concatenated
% disturbance vector as well as a bias. The resulting optimization problem is
% non-convex, and |SReachTools| formulates a difference-of-convex program to
% solve this optimization problem to a local optimum. Since affine disturbance
% feedback controllers can not satisfy hard control bounds, we relax the control
% bounds to be probabilistically violated with at most a probability of 0.01.
% After obtaining the affine feedback controller, we compute a lower bound to
% the maximal reach probability in the event saturation is applied to satisfy
% the hard control bounds. Due to its incorporation of state-feedback, this
% approach typically permits the construction of the highest underapproximative
% probability guarantee.  
if chance_affine_run_gauss
    fprintf('\n\nSReachPoint with chance-affine\n');
    opts = SReachPointOptions('term', 'chance-affine',...
        'max_input_viol_prob', 1e-2, 'verbose',2);
    timerVal = tic;
    [prob_chance_affine_gauss, opt_input_vec_chance_affine_gauss,...
        opt_input_gain_chance_affine_gauss] = SReachPoint('term', 'chance-affine',...
            sys_gauss, init_state_chance_affine_gauss, target_tube, opts);
    elapsed_time_chance_affine_gauss = toc(timerVal);
    fprintf('Computation time: %1.3f\n', elapsed_time_chance_affine_gauss);
    if prob_chance_affine_gauss > 0
        % mean_X = Z * x_0 + H * (M \mu_W + d) + G * \mu_W
        opt_mean_X_chance_affine_gauss = Z * init_state_chance_affine_gauss +...
            H * opt_input_vec_chance_affine_gauss + ...
            (H * opt_input_gain_chance_affine_gauss + G) * muW_gauss;
        % Optimal mean trajectory construction
        opt_mean_traj_chance_affine_gauss = reshape(opt_mean_X_chance_affine_gauss, ...
            sys_gauss.state_dim,[]);
        % Check via Monte-Carlo simulation
        concat_state_realization_cca_gauss = generateMonteCarloSims(...
            n_mcarlo_sims, sys_gauss, init_state_chance_affine_gauss, ...
            time_horizon, opt_input_vec_chance_affine_gauss,...
            opt_input_gain_chance_affine_gauss, 1);
        mcarlo_result = target_tube.contains(concat_state_realization_cca_gauss);
        simulated_prob_chance_affine_gauss = sum(mcarlo_result)/n_mcarlo_sims;
    else
        simulated_prob_chance_affine_gauss = NaN;
    end
    fprintf('SReachPoint underapprox. prob: %1.2f | Simulated prob: %1.2f\n',...
        prob_chance_affine_gauss, simulated_prob_chance_affine_gauss);
    fprintf('Computation time: %1.3f\n', elapsed_time_chance_affine_gauss);
end

%% Summary of results
% For ease of comparison, we list the probability estimates, the
% Monte-Carlo simulation validations, and the computation times once again.
% We also plot the mean trajectories.
figure(101);
clf;
hold on;
for itt = 0:time_horizon
    if itt==0
        % Remember the first the tube
        h_target_tube = plot(target_tube_cell{1},'alpha',0.5,'color','y');
    else
        plot(target_tube_cell{itt+1},'alpha',0.08,'LineStyle',':','color','y');
    end            
end
axis equal        
h_nominal_traj = scatter(center_box(1,:), center_box(2,:), 50,'ks','filled');        
h_vec = [h_target_tube, h_nominal_traj];
legend_cell = {'Target tube', 'Nominal trajectory'};
% Plot the optimal mean trajectory from the vertex under study
if chance_open_run_gauss
    h_opt_mean_ccc_gauss = scatter(...
          [init_state_chance_open_gauss(1), opt_mean_traj_chance_open_gauss(1,:)], ...
          [init_state_chance_open_gauss(2), opt_mean_traj_chance_open_gauss(2,:)], ...
          30, 'bo', 'filled','DisplayName', 'Mean trajectory (chance-open)');
    legend_cell{end+1} = 'Mean trajectory (chance-open)';       
    h_vec(end+1) = h_opt_mean_ccc_gauss;
    disp('>>> SReachPoint with chance-open')
    fprintf('SReachPoint underapprox. prob: %1.2f | Simulated prob: %1.2f\n',...
        prob_chance_open_gauss, simulated_prob_chance_open_gauss);
    fprintf('Computation time: %1.3f\n', elapsed_time_chance_open_gauss);    
end
if genzps_open_run_gauss
    h_opt_mean_genzps_gauss = scatter(...
          [init_state_genzps_open_gauss(1), opt_mean_traj_genzps_open_gauss(1,:)], ...
          [init_state_genzps_open_gauss(2), opt_mean_traj_genzps_open_gauss(2,:)], ...
          30, 'kd','DisplayName', 'Mean trajectory (genzps-open)');
    legend_cell{end+1} = 'Mean trajectory (genzps-open)';  
    h_vec(end+1) = h_opt_mean_genzps_gauss;
    disp('>>> SReachPoint with genzps-open')
    fprintf('SReachPoint underapprox. prob: %1.2f | Simulated prob: %1.2f\n',...
        prob_genzps_open_gauss, simulated_prob_genzps_open_gauss);
    fprintf('Computation time: %1.3f\n', elapsed_time_genzps_gauss);    
end
if particle_open_run_gauss
    h_opt_mean_particle_gauss = scatter(...
          [init_state_particle_open_gauss(1), opt_mean_traj_particle_open_gauss(1,:)], ...
          [init_state_particle_open_gauss(2), opt_mean_traj_particle_open_gauss(2,:)], ...
          30, 'r^', 'filled','DisplayName', 'Mean trajectory (particle-open)');  
    legend_cell{end+1} = 'Mean trajectory (particle-open)';    
    h_vec(end+1) = h_opt_mean_particle_gauss;
    disp('>>> SReachPoint with particle-open')
    fprintf('SReachPoint approx. prob: %1.2f | Simulated prob: %1.2f\n',...
        prob_particle_open_gauss, simulated_prob_particle_open_gauss);
    fprintf('Computation time: %1.3f\n', elapsed_time_particle_gauss);
end
if voronoi_open_run_gauss
    h_opt_mean_voronoi_open_gauss = scatter(...
          [init_state_voronoi_open_gauss(1), opt_mean_traj_voronoi_open_gauss(1,:)], ...
          [init_state_voronoi_open_gauss(2), opt_mean_traj_voronoi_open_gauss(2,:)], ...
          30, 'cv', 'filled','DisplayName', 'Mean trajectory (voronoi-open)');  
    legend_cell{end+1} = 'Mean trajectory (voronoi-open)';    
    h_vec(end+1) = h_opt_mean_voronoi_open_gauss;
    disp('>>> SReachPoint with voronoi-open')
    fprintf('SReachPoint approx. prob: %1.2f | Simulated prob: %1.2f\n',...
        prob_voronoi_open_gauss, simulated_prob_voronoi_open_gauss);
    fprintf('Computation time: %1.3f\n', elapsed_time_voronoi_gauss);    
end
if chance_affine_run_gauss
    h_opt_mean_chance_affine_gauss = scatter(...
          [init_state_chance_affine_gauss(1), opt_mean_traj_chance_affine_gauss(1,:)], ...
          [init_state_chance_affine_gauss(2), opt_mean_traj_chance_affine_gauss(2,:)], ...
          30, 'ms', 'filled','DisplayName', 'Mean trajectory (chance-affine)');
    legend_cell{end+1} = 'Mean trajectory (chance-affine)';
    h_vec(end+1) = h_opt_mean_chance_affine_gauss;
    polytopesFromMonteCarloSims(...
            concat_state_realization_cca_gauss(sys_gauss.state_dim+1:end,:), sys_gauss.state_dim,...
            [1,2], {'color','m','edgecolor','m','linewidth',2,'alpha',0.15,'LineStyle',':'});
    disp('>>> SReachPoint with chance-affine')
    fprintf('SReachPoint underapprox. prob: %1.2f | Simulated prob: %1.2f\n',...
        prob_chance_affine_gauss, simulated_prob_chance_affine_gauss);
    fprintf('Computation time: %1.3f\n', elapsed_time_chance_affine_gauss);    
end
legend(h_vec, legend_cell, 'Location','EastOutside', 'interpreter','latex');
xlabel('x');
ylabel('y');
axis equal
box on;
drawnow;

%% With beta distribution (non-Gaussian case)
n_mcarlo_sims = 1e5;                        % Monte-Carlo simulation particles
n_mcarlo_sims_affine = 1e3;                 % For affine controllers

eta_dist_nongauss = RandomVector('UserDefined', @(N) ...
    0.5*(betarnd(2, 2, 2, N) - 0.5));

[sys_nongauss, heading_vec] = getDubinsCarLtv('add-dist', turning_rate, init_heading, ...
    sampling_time, input_space, dist_matrix, eta_dist_nongauss);

% Compute the mean trajectory of the concatenated disturbance vector
muW_nongauss = sys_nongauss.dist.concat(time_horizon).mean();


% Specifying initial states and which options to run
particle_open_run_nongauss = 0;
voronoi_open_run_nongauss = 0;
voronoi_affine_run_nongauss = 1;

% Initial states for each of the method
init_state_particle_open_nongauss = [2;1];
init_state_voronoi_open_nongauss = init_state_particle_open_nongauss; %[2;2] + [-2;+0.5];
init_state_voronoi_affine_nongauss = init_state_particle_open_nongauss; %[2;2] + [-1.5;1];

%% |SReachPoint|: |particle-open| (non-Gaussian case)
% This method is discussed in <http://doi.org/10.1109/CDC.2013.6760626
% Lesser et. al., Conference on Decision and Control, 2013>.
%
% This approach implements the particle control approach to compute an open-loop
% controller. It is a sampling-based technique and hence the resulting
% probability estimate is random with its variance going to zero as the number
% of samples considered goes to infinity. Note that since a mixed-integer linear
% program is solved underneath with the number of binary variables corresponding
% to the number of particles, using too many particles can cause an exponential
% increase in computational time.
if particle_open_run_nongauss
    fprintf('\n\nSReachPoint with particle-open\n');
    opts = SReachPointOptions('term','particle-open','verbose',1,...
        'n_particles',50);
    timerVal = tic;
    [prob_particle_open_nongauss, opt_input_vec_particle_open_nongauss] =...
        SReachPoint('term', 'particle-open', sys_nongauss,...
            init_state_particle_open_nongauss, target_tube, opts);
    elapsed_time_particle_nongauss = toc(timerVal);
    if prob_particle_open_nongauss > 0
        % Optimal mean trajectory construction
        % mean_X = Z * x_0 + H * U + G * \mu_W
        opt_mean_X_particle_open_nongauss =  Z * init_state_particle_open_nongauss + ...
            H * opt_input_vec_particle_open_nongauss + G * muW_nongauss;
        opt_mean_traj_particle_open_nongauss =...
            reshape(opt_mean_X_particle_open_nongauss, sys_nongauss.state_dim,[]);
        % Check via Monte-Carlo simulation
        concat_state_realization_pao = generateMonteCarloSims(n_mcarlo_sims, ...
            sys_nongauss, init_state_particle_open_nongauss, time_horizon,...
            opt_input_vec_particle_open_nongauss);
        mcarlo_result = target_tube.contains(concat_state_realization_pao);
        simulated_prob_particle_open_nongauss = sum(mcarlo_result)/n_mcarlo_sims;
    else
        simulated_prob_particle_open_nongauss = NaN;
    end
    fprintf('SReachPoint approx. prob: %1.2f | Simulated prob: %1.2f\n',...
        prob_particle_open_nongauss, simulated_prob_particle_open_nongauss);
    fprintf('Computation time: %1.3f\n', elapsed_time_particle_nongauss);
end


%% |SReachPoint|: |voronoi-open| (non-Gaussian case)
% This method is discussed in <arxiv_link_TODO Sartipizadeh et. al., 
% American Control Conference, 2019 (submitted)>
%
% This approach implements the undersampled particle control approach to compute
% an open-loop controller. It computes, using k-means, a representative sample
% realization of the disturbance which is significantly smaller. This
% drastically improves the computational efficiency of the particle control
% approach. Further, because it uses Hoeffding's inequality, the user can
% specify an upper-bound on the overapproximation error. The undersampled
% probability estimate is used to create a lower bound of the solution
% corresponding to the original particle control problem with appropriate
% (typically large) number of particles. Thus, this has all the benefits of the
% |particle-open| option, with additional benefits of being able to specify a
% maximum overapproximation error as well being computationally tractable.
if voronoi_open_run_nongauss
    fprintf('\n\nSReachPoint with voronoi-open\n');
    opts = SReachPointOptions('term','voronoi-open','verbose',1,...
        'max_overapprox_err', 1e-2, 'n_kmeans', 50);
    timerVal = tic;
    [prob_voronoi_open_nongauss, opt_input_vec_voronoi_open_nongauss] = ...
        SReachPoint('term', 'voronoi-open', sys_nongauss,...
        init_state_voronoi_open_nongauss, target_tube, opts);
    elapsed_time_voronoi_nongauss = toc(timerVal);
    if prob_voronoi_open_nongauss > 0
        % Optimal mean trajectory construction
        % mean_X = Z * x_0 + H * U + G * \mu_W
        opt_mean_X_voronoi_open_nongauss =  Z * init_state_voronoi_open_nongauss + ...
            H * opt_input_vec_voronoi_open_nongauss + G * muW_nongauss;
        opt_mean_traj_voronoi_open_nongauss = ...
            reshape(opt_mean_X_voronoi_open_nongauss, sys_nongauss.state_dim,[]);
        % Check via Monte-Carlo simulation
        concat_state_realization_voo = generateMonteCarloSims(n_mcarlo_sims, ...
            sys_nongauss, init_state_voronoi_open_nongauss,time_horizon,...
            opt_input_vec_voronoi_open_nongauss);
        mcarlo_result = target_tube.contains(concat_state_realization_voo);
        simulated_prob_voronoi_open_nongauss = sum(mcarlo_result)/n_mcarlo_sims;
    else
        simulated_prob_voronoi_open_nongauss = NaN;
    end
    fprintf('SReachPoint approx. prob: %1.2f | Simulated prob: %1.2f\n',...
        prob_voronoi_open_nongauss, simulated_prob_voronoi_open_nongauss);
    fprintf('Computation time: %1.3f\n', elapsed_time_voronoi_nongauss);
end

%% |SReachPoint|: |voronoi-affine|
% This method extends our previous work in <https://arxiv.org/abs/1811.03643 
% Sartipizadeh, et. al., American Control Conference, 2019 (submitted)> to
% compute an affine controller. This work will be made available online
% soon. TODO
%
if voronoi_affine_run_nongauss
    fprintf('\n\nSReachPoint with voronoi-affine\n');
    opts = SReachPointOptions('term', 'voronoi-affine',...
        'max_input_viol_prob', 2e-1, 'verbose', 2, 'n_kmeans', 30,...
        'max_overapprox_err', 6e-2, 'failure_risk', 1e-4, 'bigM', 1e2);
    timerVal = tic;
    [prob_voronoi_affine_nongauss, opt_input_vec_voronoi_affine_nongauss,...
        opt_input_gain_voronoi_affine_nongauss, kmeans_info_affine_nongauss] = SReachPoint( ...
            'term', 'voronoi-affine', sys_nongauss, init_state_voronoi_affine_nongauss,...
                target_tube, opts);
    elapsed_time_voronoi_affine_nongauss = toc(timerVal);
    fprintf('Computation time: %1.3f\n', elapsed_time_voronoi_affine_nongauss);
    if prob_voronoi_affine_nongauss > 0
        % mean_X = Z * x_0 + H * (M \mu_W + d) + G * \mu_W
        opt_mean_X_voronoi_affine_nongauss = Z * init_state_voronoi_affine_nongauss +...
            H * opt_input_vec_voronoi_affine_nongauss + ...
            (H * opt_input_gain_voronoi_affine_nongauss + G) * muW_nongauss;
        % Optimal mean trajectory construction
        opt_mean_traj_voronoi_affine_nongauss = reshape(opt_mean_X_voronoi_affine_nongauss, ...
            sys_nongauss.state_dim,[]);
        % Check via Monte-Carlo simulation
        concat_state_realization_voa = generateMonteCarloSims(...
            n_mcarlo_sims_affine, sys_nongauss, init_state_voronoi_affine_nongauss,...
             time_horizon, opt_input_vec_voronoi_affine_nongauss,...
             opt_input_gain_voronoi_affine_nongauss, 1);
        mcarlo_result = target_tube.contains(concat_state_realization_voa);
        simulated_prob_voronoi_affine_nongauss = sum(mcarlo_result)/n_mcarlo_sims_affine;
    else
        simulated_prob_voronoi_affine_nongauss = NaN;
    end
    fprintf('SReachPoint approx. prob: %1.2f | Simulated prob: %1.2f\n',...
        prob_voronoi_affine_nongauss, simulated_prob_voronoi_affine_nongauss);
end



%% Summary of results
% For ease of comparison, we list the probability estimates, the
% Monte-Carlo simulation validations, and the computation times once again.
% We also plot the mean trajectories.
figure(102);
clf;
hold on;
for itt = time_horizon:-1:0
    if itt==0
        % Remember the first the tube
        h_target_tube = plot(target_tube_cell{1},'LineStyle','--','alpha',0.2,'color','y');
    else
        plot(target_tube_cell{itt+1},'alpha',0.2,'LineStyle','--','color','y');        
    end            
end
axis equal        
%h_nominal_traj = scatter(center_box(1,:), center_box(2,:), 50,'ks','filled');        
h_vec = [h_target_tube];%, h_nominal_traj];
legend_cell = {'Target tube'};%, 'Nominal trajectory'};
if particle_open_run_nongauss
    if prob_particle_open_nongauss > 0
        h_opt_mean_particle_nongauss = scatter(...
              [init_state_particle_open_nongauss(1),...
                opt_mean_traj_particle_open_nongauss(1,:)], ...
              [init_state_particle_open_nongauss(2),...
                opt_mean_traj_particle_open_nongauss(2,:)], ...
              30, 'r^', 'filled','MarkerEdgeColor','k', 'DisplayName', 'Mean trajectory (particle-open)');  
        legend_cell{end+1} = 'Mean trajectory (particle-open)';    
        h_vec(end+1) = h_opt_mean_particle_nongauss;        
        polytopesFromMonteCarloSims(...
            concat_state_realization_pao(sys_nongauss.state_dim+1:end,:), sys_nongauss.state_dim,...
            [1,2], {'color','r','edgecolor','r','linewidth',2,'alpha',0.15,'LineStyle',':'});
    end
    disp('>>> SReachPoint with particle-open')
    fprintf('SReachPoint approx. prob: %1.2f | Simulated prob: %1.2f\n',...
        prob_particle_open_nongauss, simulated_prob_particle_open_nongauss);
    fprintf('Computation time: %1.3f\n', elapsed_time_particle_nongauss);
end
if voronoi_open_run_nongauss 
    if prob_voronoi_open_nongauss > 0
        h_opt_mean_voronoi_nongauss = scatter(...
              [init_state_voronoi_open_nongauss(1), ...
                opt_mean_traj_voronoi_open_nongauss(1,:)], ...
              [init_state_voronoi_open_nongauss(2), ...
                opt_mean_traj_voronoi_open_nongauss(2,:)], ...
              30, 'cv', 'filled','MarkerEdgeColor','k', 'DisplayName', 'Mean trajectory (voronoi-open)');  
        legend_cell{end+1} = 'Mean trajectory (voronoi-open)';    
        h_vec(end+1) = h_opt_mean_voronoi_nongauss;
        polytopesFromMonteCarloSims(...
            concat_state_realization_voo(sys_nongauss.state_dim+1:end,:), sys_nongauss.state_dim,...
            [1,2], {'color','c','edgecolor','c','linewidth',2,'alpha',0.15,'LineStyle',':'});
    end
    disp('>>> SReachPoint with voronoi-open')
    fprintf('SReachPoint approx. prob: %1.2f | Simulated prob: %1.2f\n',...
        prob_voronoi_open_nongauss, simulated_prob_voronoi_open_nongauss);
    fprintf('Computation time: %1.3f\n', elapsed_time_voronoi_nongauss);    
end
if voronoi_affine_run_nongauss
    if prob_voronoi_affine_nongauss > 0
        h_opt_mean_voronoi_affine_nongauss = scatter(...
              [init_state_voronoi_affine_nongauss(1), opt_mean_traj_voronoi_affine_nongauss(1,:)], ...
              [init_state_voronoi_affine_nongauss(2), opt_mean_traj_voronoi_affine_nongauss(2,:)], ...
              30, 'bs', 'filled','DisplayName', 'Mean trajectory (voronoi-affine)');
        legend_cell{end+1} = 'Mean trajectory (voronoi-affine)';
        h_vec(end+1) = h_opt_mean_voronoi_affine_nongauss;
        polytopesFromMonteCarloSims(...
            concat_state_realization_voa(sys_nongauss.state_dim+1:end,:), sys_nongauss.state_dim,...
            [1,2], {'color','b','edgecolor','b','linewidth',2,'alpha',0.15,'LineStyle',':'});
    end
    disp('>>> SReachPoint with voronoi-affine')
    fprintf('SReachPoint approx. prob: %1.2f | Simulated prob: %1.2f\n',...
        prob_voronoi_affine_nongauss, simulated_prob_voronoi_affine_nongauss);
    fprintf('Computation time: %1.3f\n', elapsed_time_voronoi_affine_nongauss);    
end

legend(h_vec, legend_cell, 'Location','EastOutside', 'interpreter','latex');
xlabel('x');
ylabel('y');
axis equal
box on;
