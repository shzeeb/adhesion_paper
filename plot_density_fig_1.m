%plot density ABM vs PDEs for the SINGLE SPECIES SYSTEM (FIGURE 1 IN PAPER)

clear;

%set IC=3 for single species
IC=3;

%final time (from simulation)
T_final=1000;

%adhesion strenght
q=0;

%path
path=string(pwd);

%name of .mat file to read (omit extension)
full_path_simul="adhesion_IC="+num2str(IC)+"_q="+num2str(q)+"_T="+num2str(T_final);

%load data from the from named above
data_simul=load(full_path_simul+".mat");

%load determinisitc data

file_name_det="ts_2d_pde_data_q_"+num2str(q)+"_T_"+num2str(T_final)+".mat";
data_det=load(file_name_det);

%extract data to plot
rec_mat_full=data_simul.rec_mat_full;

%at which simulation time to plot ABM solution
plot_times=[0 100,1000];

%index of the above times in the solution matrix
plot_ind=[1 2 11];

%%plotting code
%define xspan vector to plot over
x=linspace(1,400,data_simul.ncols);


%go through all the plot times
for i=1:length(plot_times)

    
    %calculate averages
    % for speices 1
    rec_mat=rec_mat_full(:,:,plot_ind(i),:);
    ave_dens=mean(mean(rec_mat==1,4),1);

    
    %plot ABM solution
    plot(x,ave_dens,'Color',[0/255 0/255 0/255],'LineWidth',3);

    hold on
   
    plot(data_det.x,data_det.sol_1D(i,:),'Color',[255/255 200/255 0/255],'LineStyle','--','LineWidth',3);
    
    hold on

end

%set plot pars
xlim([120 280])
xticks([1 120 280]);
yticks([0 1]);
xlabel('x');
ylabel('density');
ax=gca;
ax.FontSize=30; 

%name figure to export
fig_name="adhesion_single_species_q_"+num2str(q);

%export figure as pdf
% exportgraphics(ax,fig_name+'.pdf')
