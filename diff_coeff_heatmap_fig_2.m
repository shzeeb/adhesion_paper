% plot diffusion coefficient heatmap ofr Figure 2
clear
close all

%diffusion coefficient
D1=@(q,u) (1-q).^(4*u).*(1-4.*u.*(u-1).*log(1-q));

u=0:0.001:1;
q=0:0.001:1;

%degine u,q mesh
[U,Q]=meshgrid(u,q);

%evaluate the diffusion coeffcient at the mesh points
D=D1(Q,U);

figure

imagesc(D,'XData',u,'YData',q)

%Split the colourmap up into 1000000 diferent gradations
cmap=parula(1000000);

%Flip it upside down
cmap=flipud(cmap);
% Get rid of the first Q entries of the colormap
Q=0;
cmap(1:Q,:)=[];
%Change the colour of the first P entries of the clour map (Where P is chosen to be the boundary os positive and negative by trial and error (i.e. we will get a different colour below zero now)) 
P=round(103731*(1000000-Q)/1000000);

cmap(1:P,:)=ones(P,1)*cmap(90000,:); %Gives a nice yellow
colormap(cmap);
colorbar

set(gca,'YDir','normal')

xlabel('C');
ylabel('q'); 
set(gca,'FontSize',20)

hold on

%plot the threshold curve above which D<0
plot(u,plot(u,1-exp(1./(4*u.*(u-1))),'k--','LineWidth',2))
