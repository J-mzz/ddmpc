%% Example 2 of "Robust Data-Driven Receding Horizon Control" (Fig. 2)
%
% Same system and constraints as Example 1 (example1.m), with a larger noise
% bound. The simulation applies a constant disturbance w_k = w (see
% utils/updateDynConstNoise.m), which is equivalent to a shifted system
% without process noise.
%
% Unknown system:  x_{k+1} = A*x_k + B*u_k + w_k,  ||w_k||_inf <= eps
% Constraints:     X = {x : ||x||_inf <= 1},  U = {u : |u| <= 1}
%
% Steps
%   0. Collect T noisy samples and build the consistency set of (A, B).
%   1. Compute the largest robust invariant set inside X (Algorithm 2 in
%      the paper), represented as W*x <= d.
%   2. Check that this set is invariant: the data-driven controller must
%      reach a contraction rate lambda <= 1 at every vertex.
%   3. Simulate from x0 without receding horizon (the consistency set is
%      fixed) and with receding horizon (the consistency set is refined
%      with every new measurement).
%
% Figures: largest invariant set, trajectories, control inputs, and
% contraction rates, with and without receding horizon (RH).

close all;clear all;clc
yalmip('clear')

%% Settings
eps    = 1e-1;      % bound on the process noise, ||w||_inf <= eps
x0     = [-1; .8];  % initial state of the closed-loop simulation
n_iter = 20;        % closed-loop steps
n_inv  = 20;        % iterations of Algorithm 2 (largest invariant set)

rng(100)
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'utils')))

%% Step 0: constraints, system and data
% state constraint X = {x : W*x <= d}
W = [1,  0;
     -1, 0;
     0,  1;
     0, -1];
d = 1*[1; 1; 1; 1];

% vertices of X are the rows of V_X
V_X = con2vert(W,d);

% input constraint U = [-u, u]
u = 1;

% ground-truth system (used to generate the data and to simulate), see utils/getSystem.m:
% A = [0, -0.99; 0.99, 0], B = [0; 1]
eg = 1;

% problem dimensions
Cons = setCons(eg);
Cons.u = u;

% noise set {w : ||V*w||_inf <= 1}
V = 1/eps*eye(2);

% T noisy samples with x in [-2, 2]^2 and u in [-1, 1]
Data = getDataRND(eg,eps,Cons.T);

%% Consistency set of (A, B) compatible with the data: N*theta <= e
[A,B,xi] = getCons(Data,Cons);   % data matrices (not the system matrices)

N = [A,  B;
    -A, -B];

% check if enough data are collected
if rank(N) ~= size(N,2)
    error('Not enough data are collected!')
end

e = [eps*ones(size(xi,1),1) + xi;
     eps*ones(size(xi,1),1) - xi];

% remove redundant constraints
[N, e] = nontrivial_constraints(N, e);

% vertices of the consistency set, one row per vertex (A_i, B_i)
V_AB = con2vert(N,e);

% append the noise set to the constraints
N = [N, zeros(size(N,1), 2);
    zeros(2,size(N,2)), V;
    zeros(2,size(N,2)), -V];

e = [e;
     ones(2,1);
     ones(2,1)];

Cons.N = N;
Cons.e = e;

%% Step 1: largest robust invariant set inside X (Algorithm 2)
figure()
polyX = polyshape(V_X(:,1), V_X(:,2));
plot(polyX)
title('Max Forward Invariant Set', 'FontSize', 12)
hold on
xlim([-1,1])
ylim([-1,1])
axis square

warning off
for ii = 1:n_inv
    % one-step robust backward reachable set R of the current set
    [W, d] = maxis(V_AB, W, d, u, eps);

    % vertices of R
    V_R = con2vert(W,d);
    polyR = polyshape(V_R(:,1), V_R(:,2));

    % update the set to R intersected with X
    polyInters = intersect(polyX,polyR);
    V_R = polyInters.Vertices;
    [W,d] = vert2con(V_R); % normalized, W*x <= 1

    plot(polyInters)
end
warning on

%% Step 2: check that the set W*x <= d is invariant
options = sdpsettings('solver', 'mosek', 'verbose', 0);

% contraction rate and input at each vertex
L = []; U = [];

for jj = 1:size(V_R,1)
    xk = V_R(jj,:)';

    out = ddmpc(xk, W, Cons, options);
    if out.sol.problem
        error('solving failed')
    end

    U = [U, out.uk];
    L = [L, out.lk];
end

if max(L) > 1 + 1e-6
    error('The set (W,d) is NOT forward invariant')
else
    fprintf('The set (W,d) IS forward invariant\n')
end

%% Step 3a: data-driven control without RH (consistency set is fixed)

% contraction rate and input log
L0 = []; U0 = [];

Traj0 = x0;

for iter = 1:n_iter
    % current state
    xk = Traj0(:,end);

    % solve for the input uk and the contraction rate lk
    out = ddmpc(xk,W,Cons,options);
    if out.sol.problem
        error('solving failed')
    end

    uk = out.uk;
    lk = out.lk;

    % apply uk under the constant disturbance; the new sample is not
    % added to the data
    xk_next = updateDynConstNoise(eps, [], xk, uk);
    xk = xk_next([5,6])';

    Traj0 = [Traj0, xk];
    U0 = [U0, uk];
    L0 = [L0, lk];
end

%% Step 3b: data-driven receding horizon control (consistency set is refined)

% contraction rate and input log
L = []; U = [];

if max(W*x0 - d) > 1e-6
    warning('The initial condition is not inside the invariant set')
end

Traj = x0;

for iter = 1:n_iter
    % current state
    xk = Traj(:,end);

    % solve for the input uk and the contraction rate lk
    out = ddmpc(xk,W,Cons,options);
    if out.sol.problem
        error('solving failed')
    end

    uk = out.uk;
    lk = out.lk;

    % apply uk under the constant disturbance and add the new sample to
    % the data
    Data = updateDynConstNoise(eps, Data, xk, uk);

    % rebuild the consistency set with all data
    [A,B,xi] = getCons(Data,Cons);
    N = [A B; -A -B];
    e = [eps*ones(size(xi,1),1)+xi; eps*ones(size(xi,1),1)-xi];
    [N, e] = nontrivial_constraints(N, e);

    % append the noise set to the constraints
    N = [N, zeros(size(N,1), 2);
        zeros(2,size(N,2)), V;
        zeros(2,size(N,2)), -V];

    e = [e;
         ones(2,1);
         ones(2,1)];

    Cons.N = N;
    Cons.e = e;

    % next state
    xk = Data(end, [5,6])';

    Traj = [Traj, xk];
    U = [U, uk];
    L = [L, lk];
end

%% Plot trajectories, inputs, and contraction rates
figure()
plot(Traj(1,:), Traj(2,:),'o--')
hold on
plot(Traj0(1,:), Traj0(2,:),'*--')
legend('RH', 'No RH', 'FontSize', 12)
title('Trajectory', 'FontSize', 12)
xlim([-1,1])
ylim([-1,1])
axis square

figure()
plot(1:length(U), U, 'o--');
hold on
plot(1:length(U0), U0, '*--');
legend('RH', 'No RH', 'FontSize', 12)
title('Control with Bound=1', 'FontSize', 12)
axis square

figure()
plot(1:length(L), L, 'o--');
hold on
plot(1:length(L0), L0, '*--');
legend('RH', 'No RH', 'FontSize', 12)
title('Contraction with Bound=1', 'FontSize', 12)
axis square
