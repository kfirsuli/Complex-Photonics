clear all; close all;clc
%Initializing Hologram Matrices
H=1080; V=1080;%%Number of Horizontal and Vertical pixels             
x=-H/2:1:(H/2-1);y=-V/2:1:(V/2-1);
[X,Y]=meshgrid(x, y);                                                  
phi=angle(X+1i*Y);%%Azimuthal angle
n=0;%%Topological charge
nx=0;ny=0;%%Number of horizontal and vertical grooves
gx=nx/H; gy=ny/V;
SR=0.5;
w0=1;%%Beam size at the SLM in  mm 
Pixel=8;%Pixel size in microns
[turb] = Turb(H,V,SR, w0,Pixel);
Hol=mod(turb+n*phi+2*pi*(Y*gy+X*gx),2*pi);
SLM=Hol/max(Hol(:))*255;
fig=figure(1);
imagesc(SLM)
colormap gray
set(fig,'Position',[120 0 H V],'MenuBar','none','ToolBar','none','resize','off');%fullscreen SLM
set(gca,'position',[0 0 1 1],'Visible','off')
axis off;
axis image
