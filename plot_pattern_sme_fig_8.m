%plots the one-dimensional SME patterns in Figure 8.

clear;

%adhesion parameters
p=0.9;
q=0.3;
r=0;
rho=1;

path=string(pwd);

file_name_det="sme_1D_rho="+num2str(rho)+"_p="+num2str(p)+"_q="+num2str(q)+"_r="+num2str(r);

%load data
data_det=load(file_name_det+".mat");

%time to plot pattern at
plot_times=[0 1000 100000];


for i=1:length(plot_times)

    A_data=squeeze(data_det.C(find(data_det.t_vec==plot_times(i)),1,:));
    B_data=squeeze(data_det.C(find(data_det.t_vec==plot_times(i)),2,:));
    
    if p==0 && q==0 && r==0.9
        figure;
        bar(data_det.x_vec,A_data,1,'FaceColor','r','FaceAlpha',1,'EdgeColor','none')
        hold on
        bar(data_det.x_vec,B_data,1,'FaceColor',[0, 166/255, 81/255],'FaceAlpha',0.7,'EdgeColor','none')
    
    else
        
    
        figure;
        stairs(data_det.x_vec,A_data,'Color',[1 0 0],'LineWidth',2); %5001
    
        hold on
    
        stairs(data_det.x_vec,B_data,'--','Color',[0, 166/255, 81/255],'LineWidth',2);

    end

    xticks([1 data_det.L]);
    yticks([0 1]);
    xlim([1 data_det.L])
    ylim([0 1])

    % pbaspect([3 2.2 1]);
    xlabel('x');
    ylabel('density');

    ax=gca;
    ax.FontSize=30; 

    
    fig_name="one_dimensional_ts_sme_p_"+num2str(p)+"_q_"+num2str(q)+"_r_"+num2str(r)+"_rho_"+num2str(rho)+"_t_"+num2str(plot_times(i));
    % saveas(ax,fig_name+".eps",'epsc');

end


%% movie

v = VideoWriter('sme_two_species_1D.avi');
v.FrameRate = 20;
v.Quality = 100;
open(v);

set(gcf,'position',[0,0,1080,1920]);

x_det=data_det.x_vec;

%looping over differnet plotting times, plot the density plots
for i=1:size(data_det.C,1)
    clf("reset")

    if p==0 && q==0 && r==0.9
        bar(data_det.x_vec,squeeze(data_det.C(i,1,:)),1,'FaceColor','r','FaceAlpha',1,'EdgeColor','none')
        hold on
        bar(data_det.x_vec,squeeze(data_det.C(i,2,:)),1,'FaceColor',[0, 166/255, 81/255],'FaceAlpha',0.7,'EdgeColor','none')
    
    else
        stairs(data_det.x_vec,squeeze(data_det.C(i,1,:)),'Color',[1 0 0],'LineWidth',2); %5001
    
        hold on
    
        stairs(data_det.x_vec,squeeze(data_det.C(i,2,:)),'--','Color',[0, 166/255, 81/255],'LineWidth',2);

    end

    xticks([1 data_det.L]);
    yticks([0 1]);
    xlim([1 data_det.L])
    ylim([0 1])

    % pbaspect([3 2.2 1]);
    xlabel('x');
    ylabel('density');

    ax=gca;
    ax.FontSize=30; 


    frame = getframe(gcf);

    writeVideo(v,frame);


end

close(v);
