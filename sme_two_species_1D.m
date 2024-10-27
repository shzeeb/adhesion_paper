%solve the 1D adhesion PDe using the method of lines
% instead of solving the PDE using finite difference to descretise time and
% space, the method of lines does not require descritisation of time. Space
% is discretised using finite difference. This way each descretisation
% point (i,j) has its own ODE to solve. The ODEs can be solved using
% one of matlab's ODE solver.

clear;

%adhesion strengths
p=0.25;
q=0;
r=0.25;

%movement rates
r_a=1;
r_b=1;

%swapping prob
rho=0.75;

%domain length
L=200;

%discretisation pts
K=L;

%time spacing for recording solution
dt=100;

%times at which to record solution
t_vec=0:dt:1000;

%space
x_vec=linspace(0,L,K);

T=length(t_vec);

%mean and st dev of the normally distrib purturbation to initial condition
mu=0.50;
sigma=0.01;

%set initial condition. Set IC=1 for species 6, IC=2 for figure 8 (d) and
%(g). IC==3 for figure 8 (a)
IC=1;

%initialise depending on the IC
if IC==1
    %middle 20% species A at density 1, the rest species B with density 0.5,
    C_0=zeros(2,K);
    C_0(1,floor(0.4*K)+1:ceil(0.6*K))=1;
    C_0(2,1:floor(0.4*K))=0.5;
    C_0(2,ceil(0.6*K)+1:K)=0.5;

elseif IC==2
    %well-mixed IC. Species A and species B uniformly at radnom density mu
    for i=1:K
        rand1=max(min(randn*sigma+mu,1),0);
        rand2=max(min(randn*sigma+mu,1),0);

        while rand1+rand2 > 1 || rand1+rand2<0

            rand1=max(min(randn*sigma+mu,1),0);
            rand2=max(min(randn*sigma+mu,1),0);

        end

        A_0(i)=rand1;
        B_0(i)=rand2;
    end
    
    C_0=[A_0;B_0];


elseif IC==3

    %half half
    C_0=zeros(2,K);
    A_0=max(min(randn(1,K/2)*sigma+mu,1),0);
    B_0=max(min(randn(1,K/2)*sigma+mu,1),0);
    C_0(1,1:K/2)=A_0;
    C_0(2,K/2+1:end)=B_0;

end


figure;
stairs(C_0(1,:),'Color',[1 0 0],'LineWidth',2);

hold on

stairs(C_0(2,:),'--','Color',[0, 166/255, 81/255],'LineWidth',2);

xticks([1 L]);
yticks([0 1]);
xlim([1 L])
ylim([0 1])

% pbaspect([3 2.2 1]);
xlabel('x');
ylabel('density');

ax=gca;
ax.FontSize=30; 


%solve the system of ODEs encoded in the RHS_sme_two_species function
[t,C]=ode15s(@(t,C)RHS_sme_two_species_1D(C,r_a,r_b,p,q,r,rho,K),t_vec,C_0);

%reshape C into suitable form to allow extraction of solution for each
%spepcies by space and time
C=reshape(C,T,2,K);


save("sme_1D_rho="+num2str(rho)+"_p="+num2str(p)+"_q="+num2str(q)+"_r="+num2str(r)+".mat");
