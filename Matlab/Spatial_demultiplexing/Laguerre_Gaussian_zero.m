% clear all; close all; clc
function SLM = Laguerre_Gaussian_zero(p,l)
%  H=1920;
%  V=1080;
 H=1024;
 V=1024;
  y=linspace(-(V/2),(V/2)-1,V);
 x=linspace(-(H/2),(H/2)-1,H);
x=x*8e-3;%%Scales the hologram in the V direction
y=y*8e-3;%%Scales the hologran in the H direction
[X,Y]=meshgrid(x,y);
LG = build_LG(p,l);
I=sqrt(LG.*conj(LG));
ph=angle(LG);
v=abs(LG)/max(max(abs(LG)));%%LG field normalized to unity 
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
%     Hol0=F.*sin(ph+2*pi*X*gx*(p+1)+2*pi*Y*gy*(l+1));
%     Hol0=F.*sin(ph+2*pi*X*gx);
Hol0=abs(ph);
    Hol=Hol0-min(Hol0(:));
%      SLM0=mod(Hol,1.17*pi);
%      SLM0=mod(Hol,pi);
SLM0=Hol;
    SLM = SLM0;
%     figure(5)
%     imagesc(SLM)
%     SLM=SLM0/max(SLM0(:))*150;
    

%      SLM0=mod(Hol,pi);
%      SLM=SLM0/max(SLM0(:))*150;

% fig=figure(2);
% scrsz = get(0,'ScreenSize')
% %set(fig,'Position',[1666 -100 1920 1080],'MenuBar','none','ToolBar','none','resize','off');%fullscreen SLM
% set(gca,'position',[0 0 1 1],'Visible','off')
% imagesc(SLM)
% colormap(gray)
% axis off
% %set(gca,'position',[0 0 1 1],'Visible','off')
% axis equal;