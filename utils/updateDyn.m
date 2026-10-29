function X = updateDyn(noise, X, xk, uk)
% UPDATEDYN  Apply one step of the true system and append the sample to X.
%   X = updateDyn(noise, X, xk, uk) computes x_{k+1} = A*xk + B*uk + noise
%   and appends the row [uk, index, xk(1), xk(2), x_{k+1}(1), x_{k+1}(2),
%   noise(1), noise(2)] to the data matrix X (pass X = [] to only get the
%   new row).

[A, B] = getSystem(1, eye(2));   % returns the matrices A and B

xk1 = A*xk + B*uk + noise;

i = size(X,1);
curr = [uk,i+1,xk(1),xk(2),xk1(1),xk1(2),noise(1),noise(2)];
X = [X; curr];

end
