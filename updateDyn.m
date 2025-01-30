function X = updateDyn(noise, X, xk, uk)
% eps: process noise bound

[A, B] = getSystem(1, eye(2));

xk1 = A*xk + B*uk + noise;

i = size(X,1);
curr = [uk,i+1,xk(1),xk(2),xk1(1),xk1(2),noise(1),noise(2)];
X = [X; curr];

end