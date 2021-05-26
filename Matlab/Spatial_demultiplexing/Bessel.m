clear all;close all
H=1920;%%Horizontal pixels             
V=1080;%%Vertical pixels
x=-H/2:1:(H/2-1);
y=-V/2:1:(V/2-1);
x=x*8e-3;%%Scales the hologram in the V direction
y=y*8e-3;%%Scales the hologran in the H direction
[X,Y]=meshgrid(x, y);                                                  
r=sqrt((X.^2 + Y.^2));
la=633e-6;
k=2*pi/la; %wavenumber
rho=sqrt(X.^2+Y.^2);
Z=X+1i*Y;
phi=angle(Z);
ni=1.5;
n=0;
alpha=.48*pi/180;
kr=k*(ni-1)*alpha;
w0=1.5;
zr=w0/(ni-1)/alpha;
D=rho<w0;
u1=exp(-1i*(-n*phi+kr*rho));%Bessel
nx=20;
ny=20;
gy=ny/(V*8e-3);
gx=nx/(H*8e-3);
hol=mod(-n*phi+kr*rho+2*pi*(X*gx+Y*gy),2*pi);
fig=figure;
imagesc(hol)
axis image
%axis off
colormap gray
set(fig,'Position',[120 0 1920 1080],'MenuBar','none','ToolBar','none','resize','off');%fullscreen SLM
set(gca,'position',[0 0 1 1],'Visible','off')
