%Initializing Hologram Matrices
H=1920; V=1080; %%Number of Horizontal and Vertical pixels
x=-H/2:1:(H/2-1); y=-V/2:1:(V/2-1);
[X,Y]=meshgrid (x, y); %%Meshgrid
x=x*8e-3;%(in mm) Scales the hologram in the V direction
y=y*8e-3;%(in mm) Scales the hologram in the H direction
[X,Y]=meshgrid(x, y);
rho=sqrt (X.^2+Y.^2);
nx=100;
ny=100;
gy=ny/(V*8e-3);
gx=nx/(H*8e-3);
lambda=0.63e-3;
k=2*pi/lambda;
kz=2*pi*sqrt(1/lambda^2-rho.^2);%

for zp=1:500;% Propagation distance in mm
    Hol=mod(kz*zp+2*pi*(X*gx+Y*gy),2*pi);
    Hol=Hol-min(Hol(:));
    SLM=Hol/max(Hol(:))*255;
    fig=figure(1);
    %set(fig,'Position',[1920 0 1920    1080],'MenuBar','none','ToolBar','none','resize','off');
    set(gca,'position',[0 0 1 1],'Visible','off')
    imagesc(SLM)
    colormap gray
    axis off;
    drawnow;
end
