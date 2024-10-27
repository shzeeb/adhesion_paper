%Figure 10: Plot density profiles or species A and species B from the
%two-dimensional ABM and PDEs

clear;

%adhesion strengths
p=0;
q=0;
r=0;

%swapping prob
rho=1;

%final time
T_final=1000;

%load stochastic ABM data
full_path_simul="two_dimensional_ts_p_"+num2str(p)+"_q_"+num2str(q)+"_r_"+num2str(r)+"_rho_"+num2str(rho);
data_simul=load(full_path_simul+".mat");

%extract occupancy matrix
rec_mat_full=data_simul.rec_mat_full;

%load deterministic (PDE)data
data_det=load("ts_2d_pde_data_p_"+num2str(p)+"_q_"+num2str(q)+"_r_"+num2str(r)+"_rho_"+num2str(rho)+"_T_"+num2str(T_final)+".mat");

%plotting times
plot_times=[0 100 1000];

%index at which to plot ABM solution
plot_ind=[1 2 11];

%index at which to plot deterministic solution
plot_ind_det=[1 2 3];

%%plotting code
%define xspan vector to plot over
% x=linspace(1,200,data_simul.ncols);

x=0.5:1:200;


for i=1:length(plot_times)
    
    %reshape rec_mat
    rec_mat=squeeze(rec_mat_full(:,:,plot_ind(i),:));
    
    %find mean density over lots of repeats
    mean_dens_1=mean(mean(rec_mat==1,3),1);
    mean_dens_2=mean(mean(rec_mat==2,3),1);

    %plot
    plot(x,mean_dens_1,'Color',[0 0 0 0.3],'LineWidth',2);
    
    hold on
    
    plot(x,mean_dens_2,'Color',[0 0 0 0.3],'LineWidth',2);


    hold on

    plot(data_det.x,data_det.sol_1D(plot_ind_det(i),:,1),'Color','r','LineWidth',3);

    hold on

    plot(data_det.x,data_det.sol_1D(plot_ind_det(i),:,2),'Color',[0, 166/255, 81/255],'LineWidth',3);

    hold on

end

%plotting params
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

annotation('arrow',X1,Y1,'Color','r','LineWidth',2);
annotation('arrow',X2,Y2,'Color',[0, 166/255, 81/255],'LineWidth',2);

fig_name="2D_pde_p_"+num2str(p)+"_q_"+num2str(q)+"_r_"+num2str(r)+"_rho_"+num2str(rho);

%save figure
% exportgraphics(ax,fig_name+'.pdf')