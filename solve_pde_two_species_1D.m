%solve_pde.m
%solves the adhesion pde system using pdepe matlab function. Saves the
%solution as matlab workspace at the end.

%By Shahzeb Raja Noureen
%Date created: 22/02/2020
%Last updated:

clear;
close all

%global vars
%discretisation step spatial
global dx;

%max x
global x_max;

%adhesion strengths
global q;

global p;

global r;

%swapping prob
global rho;

%movement rate species A
global r_a;

%movement rate species B
global r_b;

%diffusion coefficients
Da=r_a/2*dx^2;
Db=r_b/2*dx^2;

dx=0.01; %x_max/dx should be 200 as ncols=200
x_max=2;
p=0;%0:0.25:0.75;
q=0;%0:0.25:0.75;
r=0;%0:0.25:0.75;
rho=1;
r_a=1;
r_b=1;
T_final=1000;


%xspan and tspan
x = 0:dx:x_max;

%record solution every 100th time step
t = 0:100:T_final;

%symmetry parameter
m = 0;

sol = pdepe(m,@pde_two_species_1D,@pdeic_two_species_1D,@pdebc_two_species_1D,x,t);

sol_A = sol(:,:,1);
sol_B = sol(:,:,2);

%save workspace
save("adhesion_pde_rho="+num2str(rho)+"_p="+num2str(p)+"_q="+num2str(q)+"_r="+num2str(r)+".mat");