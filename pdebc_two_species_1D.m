function [pl,ql,pr,qr] = pdebc_two_species_1D(xl,ul,xr,ur,t)

%zero flux boundaries -- for pdepe boundary condition given in the form: 
%p(x,t,u)+q(x,t)f(x,t,u,du/dx)=0. For zero-flux boundaries p=0 for q=1.

pl = [0;0];
ql = [1;1];
pr = pl;
qr = ql;
end