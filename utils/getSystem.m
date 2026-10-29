function [f,g] = getSystem(eg,vars)
% GETSYSTEM  Ground-truth system x_{k+1} = f(x_k) + g*u_k used to generate data.
%   [f, g] = getSystem(eg, vars) evaluates f at vars (symbolic, sdpvar or
%   numeric). getSystem(1, eye(2)) returns the matrices A and B.

switch eg
    case 1
        % linear system of Examples 1 and 2
        A = [0, -0.99; 0.99, 0];
        f = A*vars;
        g = [0; 1];
end

end
