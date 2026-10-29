function Cons = setCons(eg)
% SETCONS  Problem dimensions for example eg.
%   Cons.n   number of states        Cons.T   number of initial samples
%   Cons.df  degree of f             Cons.dg  degree of g

switch eg
    case 1
        Cons.n 	  = 2;	    % # of state
        Cons.T    = 10;	    % # of sample
        Cons.df   = 1;      % degree of f
        Cons.dg   = 0;      % degree of g
end

end
