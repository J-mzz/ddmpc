function out = ddmpc(xk,W,Cons,options)
% DDMPC  One step of the robust data-driven receding horizon controller.
%   out = ddmpc(xk, W, Cons, options) solves the linear program
%
%       min_{uk, lk}  lk
%       s.t.  ||W*(A*xk + B*uk + w)||_inf <= lk  for all (A, B, w) in the
%             consistency set {Cons.N*theta <= Cons.e},   |uk| <= Cons.u,
%
%   using the theorem of alternatives (multipliers Y >= 0). If W*x <= 1
%   describes an invariant set, lk <= 1 certifies that x_{k+1} stays in it.
%
%   out.uk   control input          out.lk   optimal value of lk
%   out.Y    multipliers            out.sol  YALMIP diagnostics

n = Cons.n;
N = Cons.N;
e = Cons.e;
u = Cons.u;

% initialization
lk = sdpvar(1);
uk = sdpvar(1);

% W*(A*xk + B*uk + w) written as linear in theta = (A, B, w)
P2 = W*kron(eye(n), xk');
Q2 = W*kron(eye(n), uk');

N2 = [P2, Q2, W; -P2, -Q2, -W];
e2 = lk*ones(size(N2,1),1);

Y = sdpvar(size(N2,1), size(N,1),'full');

obj = lk;
F = [
    lk >= 1e-8, ...
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
