function X = updateDyn_eg2(eps, X, xk, uk)
% eps: process noise bound
% eps = 0;
rng(2025)
noise = eps*(2*rand(2,1)-1);

[A, B] = getSystem(1, eye(2));

xk1 = A*xk + B*uk + noise;

i = size(X,1);
curr = [uk,i+1,xk(1),xk(2),xk1(1),xk1(2),noise(1),noise(2)];
X = [X; curr];

end