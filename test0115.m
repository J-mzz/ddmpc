

close all;clear all;clc
yalmip('clear')

rng(100)
addpath("utils/")

%% This file contains the lastest DDMPC algorithm 0104
% Steps: 
% 0. Generate the data and the constraints sets X, U 
% 1. Use Franco's algorithm to find the max invariant set inside X
%    --> output Xk represented as Wx<=d
% 2. Check if the output polytope is forward invariant, normalize W by
%    W = inv(diag(d))*W, and solve the DDC problem at each vertex x:
%       min lambda
%       s.t norm(W*x, inf) <= lambda
%           with other constraints.
% 3. Data-driven receding horizon control for point x0 in the Xk


%% step 0: set up the state/input constraint polytope X, U
% X is given by: W * x <= d
% U is given by: U = [-u, u]

% case 1: 2d infinity norm <= 5
W = [1,  0;
     -1, 0;
     0,  1;
     0, -1];
d = 1*[1; 1; 1; 1];

% vertices of X are rows of V_X
% V_X = sortrows(con2vert(W,d)); 
V_X = con2vert(W,d);

% vertices of U
u = 1;



%% generate system and data
% system dynamics: x_k+1 = Ax_k + Bu_k

eg = 1;
% Ground truth system
% A = [0, 0.5; -.8, -.6], B = [0; 1]

sdpvar x1 x2
vars_sdp = [x1; x2];
[f_gt,g_gt] = getSystem(eg,vars_sdp);   % sdpvar

% setup
Cons = setCons(eg);
Cons.u = u;

% noise level:
% ||x_k+1 - Ax_k - Bu_k|| <= eps
eps = 4e-2;    % measurement and process noise
% It makes more sense assuming the same noise level in measure and process 

% \mathcal{V} = {Vx<=1}
V = 1/eps*eye(2);

% generate data with x in ||x|| <= 2 and u in [-1, 1]
Data = getDataRND(eg,eps,Cons.T);

%% Data-driven MPC
% Update the consistency set with new data xk and designed control u(xk)
% Then for all systems in the consistency set, find a controller online
% s.t. norm(W(Ax+Bu), inf) <= \lambda

%% data driven setup
% consistency set: Nx <= e
[A,B,xi] = getCons(Data,Cons);

% get N
N = [A,  B;
    -A, -B];

% check if enough data are collected
if rank(N) ~= size(N,2)
    error('Not enough data are collected!')
end

% get e
e = [eps*ones(size(xi,1),1) + xi;
     eps*ones(size(xi,1),1) - xi];

% remove trivial constraints
[N, e] = nontrivial_constraints(N, e);

% get vertices of the consistency set 
% rows are vertices (Ai, Bi)
V_AB = con2vert(N,e);

% update N with noise v
N = [N, zeros(size(N,1), 2);
    zeros(2,size(N,2)), V;
    zeros(2,size(N,2)), -V];

% update e with noise v
e = [e;
     ones(2,1);
     ones(2,1)];

Cons.N = N;
Cons.e = e;




%% step 1: apply Franco's algorithm to find the max invariant set
% within the state constraint X: Wx <= d

% plot the initial polytope X
figure()
polyX = polyshape(V_X(:,1), V_X(:,2));
plot(polyX)
title('Max Forward Invariant Set', 'FontSize', 12)
hold on
xlim([-1,1])
ylim([-1,1])
axis square

% V_AB = [0, 0.5, -.8, -.6, 0, 1]; % ground truth
warning off
for ii = 1:20 % # of iterations in Franco's algorithm
    % pause(0.1)

    % calculate the R in step 3 of Franco's
    [W, d] = maxis(V_AB, W, d, u, eps);

    % get vertices of R
    V_R = con2vert(W,d);
    polyR = polyshape(V_R(:,1), V_R(:,2));
    % plot(polyR)

    % find the intersection of R and Xk
    polyInters = intersect(polyX,polyR);
    % polyX = polyInters;          % doesn't matter
    V_R = polyInters.Vertices;

    % update Xk = R intersect Xk
    [W,d] = vert2con(V_R); % normalized

    % plot the intersection polytope
    plot(polyInters)
end
warning on

% generate gif of iteration
% exportgraphics(gcf, gifFile, Append=true);
% exportgraphics(gca, gifFile, Append=true);

%% step 2: check if the polytope is forward invariant
% data-driven control, for each vertices of X
options = sdpsettings('solver', 'mosek', 'verbose', 0);

% lambda and input log
L = []; U = [];

for jj = 1:size(V_R,1)
    % for each vertex of Xk
    xk = V_R(jj,:)';

    % solve ddmpc for lk and uk
    out = ddmpc(xk, W, Cons, options);
    if out.sol.problem
        error('solving failed')
    end
    
    U = [U, out.uk];
    L = [L, out.lk];
end

if max(L) > 1 + 1e-6
    error('The set (W,d) is NOT forward invariant\n')
else
    fprintf('The set (W,d) IS forward invariant\n')
end

%% intermediate step: Data-driven control without RH
n_iter = 50; % trajectory length.
% rng(1234)
noise = eps*(2*rand(2, n_iter)-1);

% lambda and input log
L0 = []; U0 = [];

% check the initial condition inside of x0
x0 = [1; .8]; %[-1; .8];

% data-driven receding horizon
Traj0 = x0;

for iter = 1:n_iter % # of times updating consistency set
    
    % current state
    xk = Traj0(:,end);
    
    % solve ddmpc for lk and uk
    out = ddmpc(xk,W,Cons,options);
    if out.sol.problem
        error('solving failed')
    end

    % extract solution
    uk = out.uk;
    lk = out.lk;

    % update consistency set
    xk_next = updateDyn(noise(:,iter), [], xk, uk);
    
    % get the next step
    xk = xk_next([5,6])';

    % update log
    Traj0 = [Traj0, xk];
    U0 = [U0, uk];
    L0 = [L0, lk];
end


%% step 3: Data-driven receding horizon control for point x0 in the Xk

% lambda and input log
L = []; U = [];

% check the initial condition inside of x0
% x0 = [2; 5];
if max(W*x0 - d) > 1e-6
    warning('The initial condition is not inside the invariant set')
end

% data-driven receding horizon
Traj = x0;

for iter = 1:n_iter % # of times updating consistency set
    
    % current state
    xk = Traj(:,end);
    
    % solve ddmpc for lk and uk
    out = ddmpc(xk,W,Cons,options);
    if out.sol.problem
        error('solving failed')
    end

    % extract solution
    uk = out.uk;
    lk = out.lk;

    % update consistency set
    Data = updateDyn(noise(:,iter), Data, xk, uk);
    [A,B,xi] = getCons(Data,Cons);
    N = [A B; -A -B];
    e = [eps*ones(size(xi,1),1)+xi; eps*ones(size(xi,1),1)-xi];
    [N, e] = nontrivial_constraints(N, e);

    % update N with noise v
    N = [N, zeros(size(N,1), 2);
        zeros(2,size(N,2)), V;
        zeros(2,size(N,2)), -V];
    
    % update e with noise v
    e = [e;
         ones(2,1);
         ones(2,1)];
    
    Cons.N = N;
    Cons.e = e;
    
    % get the next step
    xk = Data(end, [5,6])';

    % update log
    Traj = [Traj, xk];
    U = [U, uk];
    L = [L, lk];
end

% plot trajectory, control, and lambda
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
legend('RH', 'No RH', 'FontSize',12)
title('Contraction with Bound=1', 'FontSize',12)
axis square