function X = updateDynConstNoise(eps, X, xk, uk)
% UPDATEDYNCONSTNOISE  Same as updateDyn, with a constant disturbance.
%   X = updateDynConstNoise(eps, X, xk, uk) draws the disturbance
%   w = eps*(2*rand(2,1)-1) after reseeding the random number generator,
%   so every call applies the same w (Example 2 in the paper).
%   Note: this resets the global random stream with rng(2025).

rng(2025)
noise = eps*(2*rand(2,1)-1);

[A, B] = getSystem(1, eye(2));   % returns the matrices A and B

xk1 = A*xk + B*uk + noise;

i = size(X,1);
curr = [uk,i+1,xk(1),xk(2),xk1(1),xk1(2),noise(1),noise(2)];
X = [X; curr];

end
