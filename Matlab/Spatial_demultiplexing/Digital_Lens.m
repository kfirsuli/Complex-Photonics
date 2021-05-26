clc; clear all; close all
%Initializing Hologram Matrices
H=1920;%%Horizontal pixels             
V=1080;%%Vertical pixels
x=-H/2:1:(H/2-1);
y=-V/2:1:(V/2-1);
x=x*8e-3;
y=y*8e-3;
[X,Y]=meshgrid(x, y);  
r=sqrt(X.^2 + Y.^2);
nx=100;
ny=0;
gy=ny/(V*8e-3);
gx=nx/(H*8e-3);
lambda=0.63e-3; % units mm
k=2*pi/lambda;
ff=10;  
T=pi/lambda/ff*(X.^2+Y.^2);
HOL=mod(T+2*pi*(X*gx+Y*gy),2*pi);
SLM=HOL-min(HOL(:));
SLM=SLM/max(SLM(:))*255;
fig=figure;
set(fig,'Position',[120 0 1920 1080],'MenuBar','none','ToolBar','none','resize','off');%fullscreen SLM
set(gca,'position',[0 0 1 1],'Visible','off')
imagesc(SLM)
axis image
axis off
colormap gray


