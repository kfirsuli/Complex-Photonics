%Initializing Hologram Matrices
H=1920; V=1080;%%Number of Horizontal and Vertical pixels
x=-H/2:1:(H/2-1);y=-V/2:1:(V/2-1);
[X,Y]=meshgrid(x, y);
phi=angle(X+1i*Y);%%Azimuthal angle
l=1;%%Topological charge
nx=100;ny=100;%%Number of horizontal and vertical grooves
gx=nx/H; gy=ny/V;
Hol=mod(l*phi+2*pi*(Y*gy+X*gx),2*pi);
%%%%%% Grascale normalization from [0, 2Pi]to [0 255]
SLM=Hol/max(Hol(:))*255;
fig=figure(1);
imagesc(SLM)
colormap gray
%set(fig,'Position',[1920 0 1920 1080],'MenuBar','none','ToolBar','none','resize','off');%fullscreen SLM
set(gca,'position',[0 0 1 1],'Visible','off')
axis off;