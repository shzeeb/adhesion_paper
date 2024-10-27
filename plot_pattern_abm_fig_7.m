%plotting code for Figure 7. Produces patterns from one-dimensional ABM
%data. Data can be simulated using the one-dimensional ABM code and by
%setting the appropriate IC.

clear;

%%parameters

%final time
T_final=100000;

%adhesion parameter p
p=0.9;

%adhesion parameter 1
q=0.3;

%adhesion parameter r
r=0;

%swapping probability
rho=1;

%current working directory
% path=string(pwd);

%name of .mat file to load
file_name_simul="pattern_two_species_w_prolif_1D_p="+num2str(p)+"_q="+num2str(q)+"_r="+num2str(r)+"_rho="+num2str(rho)+"_T="+num2str(T_final);%+"_200reps";

%load data
data_simul=load(file_name_simul+".mat");

%extract the recording matrix
rec_mat_full=squeeze(data_simul.rec_mat_full);

%plotting times
plot_times=data_simul.rec_times;

%adjust plotting times if plot_times(1)~=0
if plot_times(1)~=0
    plot_times=[0 plot_times];
end

%plotting indices
plot_idx=[1 11 length(plot_times)];


%looping over differnet plotting times, plot the density plots
for i=1:length(plot_idx)

    rec_mat=squeeze(rec_mat_full(:,plot_idx(i)))';

    figure;
    imagesc(rec_mat);
    pbaspect([10 1 1])
    c=zeros(3,3);

    c(1,:)=[1 1 1];
    c(2,:)=[1 0 0];
    c(3,:)=[0 166/255 81/255];

    colormap(c)

    xlim([1 data_simul.ncols]);
    xticks([1 data_simul.ncols]);
    yticks([]);
    xlabel('x');
    ax=gca;
    ax.FontSize=16;


    %name and save figure as eps
    fig_name="pattern_abm_p_"+num2str(p)+"_q_"+num2str(q)+"_r_"+num2str(r)+"_rho_"+num2str(rho)+"_t_"+num2str(plot_times(plot_idx(i)));
    % saveas(ax,fig_name+".eps",'epsc');
    
end



%% movie

v = VideoWriter('ms_adh_1D.avi');
v.FrameRate = 20;
v.Quality = 100;
open(v);

set(gcf,'position',[0,0,1080,1920]);


%looping over differnet plotting times, plot the density plots
for i=1:size(rec_mat_full,2)
    clf("reset")

    A=squeeze(rec_mat_full(:,i));

    imagesc(A')
    
    pbaspect([10 1 1])

    colormap(c)

    % ylim([0 1]);
    xlim([1 data_simul.ncols]);
    xticks([1 data_simul.ncols]);
    yticks([]);
    xlabel('x');
    % ylabel('density');
    ax=gca;
    ax.FontSize=16;


    title("t="+num2str(i*100))

    frame = getframe(gcf);

    writeVideo(v,frame);


end

close(v);
