function [Fout, gout] = maxis(V_AB, F, g, u, eps)
% MAXIS  One iteration of the largest robust invariant set computation
%   (Franco's algorithm, Algorithm 2 in the paper).
%   [Fout, gout] = maxis(V_AB, F, g, u, eps) returns a polytope
%   {x : Fout*x <= gout} of states that every system in the consistency
%   set can steer into {x : F*x <= g} with an input in [-u, u], robustly
%   to ||w||_inf <= eps (see below).
%
%   V_AB  vertices of the consistency set, one row [A(1,:), A(2,:), B']
%         per vertex
%   For each extreme input -u and u, the preimages of the set under all
%   vertices (A_i, B_i) are intersected; the result is the convex hull of
%   both intersections.

% tighten the set by the worst-case noise
d_vert = eps*[1 1 -1 -1; 1 -1 1 -1];
maxFd = max(F*d_vert, [], 2);
g = g - maxFd;

%% input u_k = -u

feas_u = -u;

for ii = 1:size(V_AB,1)
    V_xkold = [];
    % for each vertex, reshape Ai and Bi
    Ai = reshape(V_AB(ii, 1:4), 2, 2)';
    Bi = V_AB(ii, 5:6)';

    % for vertices uj in U, calculate xk s.t. F(Axk + Bu)  <= g
    % --> FA xk <= g - FB uj
    for jj = 1:length(feas_u)
        F_stack = F*Ai;
        g_stack = g-F*Bi*feas_u(jj);
        xk_old = con2vert(F_stack,g_stack);
        V_xkold = [V_xkold; xk_old];
    end

    if ii == 1
        poly1 = polyshape(V_xkold(:,1), V_xkold(:,2));
    else
        poly2 = polyshape(V_xkold(:,1), V_xkold(:,2));
        polyout = intersect(poly1,poly2);
        V = polyout.Vertices;
        V = rmmissing(V);
        if size(V, 1) == 0
            continue
        end
        [Ftemp, gtemp] = vert2con(V);
        V = con2vert(Ftemp,gtemp);
        poly1 = polyshape(V(:,1), V(:,2));
    end

end
V1 = V;

%% input u_k = u

feas_u = u;

for ii = 1:size(V_AB,1)
    V_xkold = [];
    % for each vertex, reshape Ai and Bi
    Ai = reshape(V_AB(ii, 1:4), 2, 2)';
    Bi = V_AB(ii, 5:6)';

    % for vertices uj in U, calculate xk s.t. F(Axk + Bu)  <= g
    % --> FA xk <= g - FB uj
    for jj = 1:length(feas_u)
        F_stack = F*Ai;
        g_stack = g-F*Bi*feas_u(jj);
        xk_old = con2vert(F_stack,g_stack);
        V_xkold = [V_xkold; xk_old];
    end

    if ii == 1
        poly1 = polyshape(V_xkold(:,1), V_xkold(:,2));
    else
        poly2 = polyshape(V_xkold(:,1), V_xkold(:,2));
        polyout = intersect(poly1,poly2);
        V = polyout.Vertices;
        V = rmmissing(V);
        if size(V, 1) == 0
            continue
        end
        [Ftemp, gtemp] = vert2con(V);
        V = con2vert(Ftemp,gtemp);
        poly1 = polyshape(V(:,1), V(:,2));
    end

end

V2 = V;
V = [V1; V2];

[Fout,gout] = vert2con(V);

end
