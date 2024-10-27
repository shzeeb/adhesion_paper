%plot 2D ABM engulfment and checkerboard patterns in Figure 11 of the
%paper.

%load data

p=0;
q=0;
r=0.8;
rho=1;
T_final=100000;

file_abm=load("pattern_abm_2D_p_"+num2str(p)+"_q_"+num2str(q)+"_r_"+num2str(r)+"_rho_"+num2str(rho)+".mat");

rec_mat_full=file_abm.rec_mat_full;


pattern_end=rec_mat_full(:,:,end);

cMap = interp1([0;1],[1 0 0;0 166/255 81/255],linspace(0,1,256));


t_plt=[0 10000 T_final];%86000, 50000

for i=1:length(t_plt)
    figure;
    if t_plt(i)>=100
    
        imagesc(squeeze(rec_mat_full(:,:,t_plt(i)/100)))
    
    elseif t_plt(i)==0
    
        imagesc(file_abm.A_init);
    
    end
    
    colormap(cMap)
    set(gca,'YDir','normal')
    set(gca,'FontSize',30)
    % xlim([1 file_sme.ncols]);
    % ylim([1 file_sme.nrows]);
    xticks([1 file_abm.ncols])
    yticks([1 file_abm.nrows])
    axis square
    xlabel('x');
    ylabel('y')

    % if t_plt(i)==0
    %     exportgraphics(gca,'abm_check_IC.eps') 
    % 
    % else
    %     exportgraphics(gca,"abm_check_t_"+num2str(t_plt(i))+".eps") 
    % end

end




%%
v = VideoWriter('test_2D_abm_pattern.avi');
v.FrameRate = 20;
v.Quality = 100;
open(v);

set(gcf,'position',[0,0,1080,1920]);

c=zeros(3,3);
% c(1,:)=[1 1 1];
c(1,:)=[1 0 0];
c(2,:)=[0 166/255 81/255];


for i=1:size(rec_mat_full,3)                
    imagesc(rec_mat_full(:,:,i));

    colormap(c);
    
    axis square
    pbaspect([1 1 1]);
    title_ = sprintf('t = %f',file_abm.rec_times(i));
    title(title_);
    
    frame = getframe(gcf);
    writeVideo(v,frame);

end
