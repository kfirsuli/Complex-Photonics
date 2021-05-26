clear all; close all;clc
%Initializing Hologram Matrices
H=1920; V=1080;%%Number of Horizontal and Vertical pixels             
x=-H/2:1:(H/2-1);y=-V/2:1:(V/2-1);
x=x*20e-3;%%Scales the hologram in the V direction
y=y*20e-3;%%Scales the hologran in the H direction
[X,Y]=meshgrid(x, y);                                                  
A= (X.^3+Y.^3)/3;
nx=100;ny=100;%%Number of horizontal and vertical grooves
gx=nx/H; gy=ny/V;
H=mod(A,2*pi);
%%%%%% Grayscale normalization from [0, 2Pi]to [0 255]
SLM=H/max(H(:))*255;
fig=figure(1);
imagesc(SLM)
colormap gray
axis image
%set(fig,'Position',[1920 0 1920 1080],'MenuBar','none','ToolBar','none','resize','off');%fullscreen SLM
set(gca,'position',[0 0 1 1],'Visible','off')
axis off;


