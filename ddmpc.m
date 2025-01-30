function out = ddmpc(xk,W,Cons,options)
% given consistency set and current state
% find the data-driven super-stability controller

n = Cons.n;
N = Cons.N;
e = Cons.e;
u = Cons.u;

% initialization
lk = sdpvar(1);
uk = sdpvar(1);

% paranmeter W
P2 = W*kron(eye(n), xk');
Q2 = W*kron(eye(n), uk');

N2 = [P2, Q2, W; -P2, -Q2, -W];
e2 = lk*ones(size(N2,1),1);

Y = sdpvar(size(N2,1), size(N,1),'full');

obj = lk;
F = [
    lk >= 1e-8, ... 1e-8, ...
    Y*N == N2, ...
    Y*e <= e2,
    Y >= 0, ...
    uk <= u, ...
    uk >= -u
    ];

sol = optimize(F, lk, options);

out.sol = sol;
out.lk = value(lk);
out.uk = value(uk);
out.Y = value(Y);

end