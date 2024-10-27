%plots ABM and DME densities of the one dimensional two-species model to
%create Figure 6. Loads stochastic and DME data in stored in .mat format
% and plots the densities in a single figure.

clear;

%simulation run time
T_final=1000;

%adhesion strenghts
p=0.25;
q=0;
r=0.25;

%swapping prob
rho=0.75;

%load ABM data
full_path_simul="one_dimensional_ts_p="+num2str(p)+"_q="+num2str(q)+"_r="+num2str(r)+"_rho="+num2str(rho)+"_T="+num2str(T_final);
data_simul=load(full_path_simul+".mat");

%extract occupancy matrix from ABM data
rec_mat_full=data_simul.rec_mat_full;

%load DME data
file_name_det="sme_1D_rho="+num2str(rho)+"_p="+num2str(p)+"_q="+num2str(q)+"_r="+num2str(r);
data_det=load(file_name_det+".mat");

%time index at which to plot density profile
plot_times=[0 100 1000];
plot_ind=[1 2 11];

plot_ind_det=[1 2 11];

%x vector to plot over
x=1:1:200;
x_det=linspace(1,200,200);

%reshape the occupancy matrix
rec_mat_full=squeeze(rec_mat_full);


for i=1:length(plot_times)
    
    %calculate averages
    % for speices 1
    
    rec_mat=rec_mat_full(:,plot_ind(i),:); %6
    rec_mat=squeeze(rec_mat);
    
    %calculate mean density over all the repeats contained in rec_mat
    mean_dens_1=mean(rec_mat==1,2);
    mean_dens_2=mean(rec_mat==2,2);

    %plot the density
    plot(x_det,squeeze(data_det.C(plot_ind_det(i),1,:)),'Color',[1 0 0],'LineWidth',2); %5001
    
    hold on
    
    plot(x_det,squeeze(data_det.C(plot_ind_det(i),2,:)),'Color',[0, 166/255, 81/255],'LineWidth',2);
    
    hold on
    
    plot(x,mean_dens_1,'Color',[0 0 0],'LineWidth',1);

    hold on

    plot(x,mean_dens_2,'Color',[0 0 0],'LineWidth',1);

    hold on

end


%plotting pars
xticks([1 200]);
yticks([0 1]);
xlim([1 data_simul.ncols])

% pbaspect([8 1 1]);
xlabel('x');
ylabel('density');
ax=gca;
ax.FontSize=30; 
% ax.XColor=[77/255, 77/255, 79/255];
% ax.YColor=[77/255, 77/255, 79/255];

%draws an arrow pointing in the direction of increasing time
X1 = [0.52 0.52];
Y1 = [0.92 0.70];

X2=[0.52 0.52];
Y2=[0.2 0.4];

annotation('arrow',X1,Y1,'Color','red','LineWidth',2);
annotation('arrow',X2,Y2,'Color',[0, 166/255, 81/255],'LineWidth',2);

%name fig
fig_name="one_dimensional_ts_sme_p="+num2str(p)+"_q="+num2str(q)+"_r="+num2str(r)+"_rho="+num2str(rho);

%save fig
% exportgraphics(ax,fig_name+'.pdf')
