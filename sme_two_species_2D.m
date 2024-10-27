%solve the DMEs

clear;

tic()

%adhesion strengths
p=0.2;
q=0.9;
r=0;

%movement rates
r_a=1;
r_b=1;

%swapping probability
rho=1;

%num of cols and rows
ncols=20;
nrows=20;

% size of time step
dt=100;

%vector of time to store DME solutions
t_vec=0:dt:100000;

%position vectors to evaluate solution
x_vec=linspace(0,1,ncols);
y_vec=linspace(0,1,nrows);

T=length(t_vec);

%mean of the random normally distributed perturbation to the initial state
mu=0.7;

%standard deviation of the random normally distributed perturbation to the
%initial state
sigma=0.01;
    
%initial varialbe to store initial states of domain
C_0=zeros(2,nrows,ncols);

for i=1:nrows
    for j=1:ncols
        rand1=max(min(randn*sigma+mu,1),0);

        A_0(i,j)=rand1;
        B_0(i,j)=1-rand1;
    end
end

C_0(1,:,:)=A_0;
C_0(2,:,:)=B_0;

% solve the equations using ode15s of MATLAB

[t,C]=ode15s(@(t,C)RHS_2D_sme_two_species(C,r_a,r_b,p,q,r,rho,nrows,ncols),t_vec,C_0);

%reshape the resulting vector to the correct size to extract density of
%species A and B from it
C=reshape(C,T,2,nrows,ncols);
            
% extract the solution at the final time point
% A_end=squeeze(C(end,1,:,:));
% B_end=squeeze(C(end,2,:,:));
% 
% plot it
% [x_mesh,y_mesh]=meshgrid(x_vec,y_vec);

%save variables
save("pattern_2D_sme_"+num2str(nrows)+"by"+num2str(ncols)+"_p_"+num2str(p(i))+"_q_"+num2str(q(j))+"_r_"+num2str(r(k))+"_rho_"+num2str(rho)+"_T_"+num2str(t_vec(end))+".mat");


toc()
