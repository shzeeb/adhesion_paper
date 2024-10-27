% creates the SME engulfment and checkerboard patterns shown in Figure 11
% of the paper.

clear

file_sme=load("pattern_2D_sme_20by20_p_0.2_q_0.9_r_0_rho_1_T_100000.mat");

% extract simulated data
A=squeeze(file_sme.C(:,1,:,:));
B=squeeze(file_sme.C(:,2,:,:));
%define colour map
cMap = interp1([0;1],[1 0 0;0 166/255 81/255],linspace(0,1,256));

%times to plot pattern for
t_plt=[0 10000 100000];

for i=1:length(t_plt)
    figure;
    
    
    if t_plt(i)>=100
        imagesc(squeeze(A(t_plt(i)/100,:,:)))
    
    elseif t_plt(i)==0
        imagesc(squeeze(A(1,:,:)))
    
    end  
    
    colormap(flipud(cMap))
    
    set(gca,'YDir','normal')
    set(gca,'FontSize',30)
    xticks([1 file_sme.ncols])
    yticks([1 file_sme.nrows])
    h=gca;
    h.XAxis.TickLength = [0 0];
    h.YAxis.TickLength = [0 0];
    xlabel('x');
    ylabel('y')
    clim([0 1])
    axis square

    if t_plt(i)==t_plt(end)
        cbh = colorbar ; %Create Colorbar
        cbh.Ticks = linspace(0, 1, 3) ;
        cbh.TickLabels = num2cell(["B", "A and B" "A"]);
        cbh.Ruler.TickLabelRotation=90;
    end
    
    if t_plt(i)==0
        % exportgraphics(gca,'sme_engulf_IC.eps') 

    else
        % exportgraphics(gca,"sme_engulf_t_"+num2str(t_plt(i))+".eps") 
    end
end
%% creates a movie of the two species simulated data provided

v = VideoWriter('ts_adh_2D_sme_pattern.avi');
v.FrameRate = 20;
v.Quality = 100;
open(v);

set(gcf,'position',[0,0,1080,1920]);

for i=1:size(A,1)                
    imagesc(squeeze(A(i,:,:)));

    colormap(flipud(cMap))
    
    set(gca,'YDir','normal')
    set(gca,'FontSize',30)
    xticks([1 file_sme.ncols])
    yticks([1 file_sme.nrows])
    h=gca;
    h.XAxis.TickLength = [0 0];
    h.YAxis.TickLength = [0 0];
    xlabel('x');
    ylabel('y')
    clim([0 1])
    axis square
    
    title_ = sprintf('t = %f',file_sme.t(i));
    title(title_);
    
    frame = getframe(gcf);
    writeVideo(v,frame);

end
