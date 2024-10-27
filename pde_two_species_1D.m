function [c,f,s] = pde_two_species_1D(x,t,u,dudx)

global dx;
global rho;
global p;
global q;
global r;
global Pm;
global Px;


Da=Pm/2*dx^2;
Db=Px/2*dx^2;

c = [1;1];

%flux term
f=[Da*((1-p)^(2*u(1))*(1-r)^(2*u(2))*((2*log(1-r)*(1-u(1)-u(2))+1)*u(1)*dudx(2)+...
    (2*log(1-p)*u(1)*(1-u(1)-u(2))+1-u(2))*dudx(1)))+...
    (Da+Db)*rho*(1-p)^(2*u(1))*(1-q)^(2*u(2))*(1-r)^(2*u(1)+2*u(2))*(u(2)*(1-2*log(1-r)*u(1)+2*log(1-p)*u(1))*dudx(1)-...
    u(1)*(1-2*log(1-q)*u(2)+2*log(1-r)*u(2))*dudx(2));
    
    Db*((1-q)^(2*u(2))*(1-r)^(2*u(1))*((2*log(1-r)*(1-u(1)-u(2))+1)*u(2)*dudx(1)+...
    (2*log(1-q)*u(2)*(1-u(1)-u(2))+1-u(1))*dudx(2)))+...
    (Da+Db)*rho*(1-q)^(2*u(2))*(1-q)^(2*u(1))*(1-r)^(2*u(1)+2*u(2))*(u(1)*(1-2*log(1-r)*u(2)+2*log(1-q)*u(2))*dudx(2)-...
    u(2)*(1-2*log(1-p)*u(1)+2*log(1-r)*u(1))*dudx(1))];


%no source term
s = [0;0];

end