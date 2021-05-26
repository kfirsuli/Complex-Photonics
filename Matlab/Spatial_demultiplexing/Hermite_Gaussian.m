clear all; close all; clc
 H=1920;
 V=1080;
 y=linspace(-(V/2),(V/2)-1,V);
 x=linspace(-(H/2),(H/2)-1,H);
x=x*8e-3;%%Scales the hologram in the V direction
y=y*8e-3;%%Scales the hologran in the H direction
[X,Y]=meshgrid(x,y);
phi=angle(X+1i*Y);
rho=sqrt(X.^2+Y.^2);
lambda=633e-6;
w0=1;
z=.00001;
k=2*pi/lambda;
zr=pi*w0^2/lambda;
w=w0*sqrt(1+(z/zr)^2);
R=z*(1+(zr/z)^2);
m=3;
n=4;
Hx=polyval(Hermite(m),sqrt(2)*X/w);
Hy=polyval(Hermite(n),sqrt(2)*Y/w);
rc=sqrt(2^(1-n-m)/(pi*factorial(n)*factorial(m)))/w;
HG=rc*Hx.*Hy.*exp(1i*(n+m+1)*atan(z/zr)).*exp(-rho.^2/w^2).*exp(-1i*k*rho.^2/(2*R))*exp(1i*k*z);
HG=HG/sqrt(sum(sum(abs(HG).^2)));
I=HG.*conj(HG);
ph=angle(HG);
v=abs(HG)/max(max(abs(HG)));%%LG field normalized to unity 
aux=round(v*800+1);
load fx2.mat;  
% 
  for mh=1:V;
      for nh=1:H;
          temp=aux(mh,nh);
          F(mh,nh)=fx(temp);                                
      end                                                                    
  end 
nx=50;
ny=50;
gy=ny/(V*8e-3);
gx=nx/(H*8e-3);
    Hol=F.*sin(ph+2*pi*(X*gx+Y*gy));
    Hol=Hol-min(Hol(:));
    SLM=Hol/max(Hol(:))*255;
imagesc(SLM)
colormap(gray)
fig=figure(1);
%set(fig,'Position',[190 0 1920 1080],'MenuBar','none','ToolBar','none','resize','off');%fullscreen SLM
set(gca,'position',[0 0 1 1],'Visible','off')
axis image;


