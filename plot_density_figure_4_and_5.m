%Code for Figures 4 and 5: plot the density profiles of species A and B
%from the one-dimensional ABM and the one-dimension corresponding PDEs.

clear;

% max time
T_final=1000;

%adhesion strengths
p=0;
q=0;
r=0;

%swapping prob
rho=1;

% load ABM simulated dfata
full_path_simul="one_dimensional_ts_p="+num2str(p)+"_q="+num2str(q)+"_r="+num2str(r)+"_rho="+num2str(rho)+"_T="+num2str(T_final);
data_simul=load(full_path_simul+".mat");

%extract occupancy matrix
rec_mat_full=data_simul.rec_mat_full;

%load PDE data
file_name_det="adhesion_pde_rho="+num2str(rho)+"_p="+num2str(p)+"_q="+num2str(q)+"_r="+num2str(r);
data_det=load(file_name_det+".mat");

%plotting times
plot_times=[0 100 1000];

%index of plotting time in recording matrix
plot_ind=[1 2 11];

%x span
x=linspace(0,200,200);


%reshape ocupancy amtrix
rec_mat_full=squeeze(rec_mat_full);


figure;
for i=1:length(plot_times)
    
    rec_mat=rec_mat_full(:,plot_ind(i),:);
    rec_mat=squeeze(rec_mat);

    %calculate averaged density over all the repeats
    mean_dens_A=mean(rec_mat==1,2);
    mean_dens_B=mean(rec_mat==2,2);

    %plot determistic density
    plot(data_det.x*100,data_det.sol_A(plot_ind(i),:),'Color',[1 0 0],'LineWidth',3); %5001

    hold on

    plot(data_det.x*100,data_det.sol_B(plot_ind(i),:),'Color',[0, 166/255, 81/255],'LineWidth',3);

    hold on
    
    %plot stochastic density
    plot(x,mean_dens_A,'Color',[0 0 0],'LineWidth',2);

    hold on
    
    plot(x,mean_dens_B,'Color',[0 0 0],'LineWidth',2);

    hold on

end

xticks([1 200]);
yticks([0 1]);
xlabel('x');
ylabel('density');
ax=gca;
ax.FontSize=30; 


%draw an arrow pointing downward (showing direction of time)
X1 = [0.52 0.52];
Y1 = [0.92 0.70];

X2=[0.52 0.52];
Y2=[0.2 0.4];

annotation('arrow',X1,Y1,'Color','red','LineWidth',2);
annotation('arrow',X2,Y2,'Color',[0, 166/255, 81/255],'LineWidth',2);



fig_name="one_dimensional_ts_p="+num2str(p)+"_q="+num2str(q)+"_r="+num2str(r)+"_rho="+num2str(rho);

% exportgraphics(ax,fig_name+'.pdf')
