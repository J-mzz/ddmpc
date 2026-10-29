function x = getDataRND(eg,eps,N)
% GETDATARND  N noisy samples of the true system at random states.
%   x = getDataRND(eg, eps, N) samples x in [-2, 2]^2 and u in [-1, 1] and
%   returns one row [u, i, x1, x2, x1+, x2+, w1, w2] per sample, where
%   x+ = f(x) + g*u + w and ||w||_inf <= eps (see getSystem.m).

box = 2;

syms x1 x2
vars = [x1;x2];
[f,g] = getSystem(eg,vars);

x = [];
for i = 1:N
    
    u = 2*rand(1,1)-1;

	noise = eps*(2*rand(2,1)-1);
    
    x0 = (2*rand(2,1)-1) * box;
    
    x1 = x0(1);
    x2 = x0(2);
    
    xdot = subs(f+g*u+noise);
    
    curr = [u,i,x1,x2,xdot(1),xdot(2), noise(1), noise(2)];
    
    x = [x;curr];
end

x = double(x);
end