% clear all; close all; clc
function SLM = Linearly_polarized(i)
 H=1024;
 V=1024;
  y=linspace(-(V/2),(V/2)-1,V);
 x=linspace(-(H/2),(H/2)-1,H);
x=x*8e-3;%%Scales the hologram in the V direction
y=y*8e-3;%%Scales the hologran in the H direction
[X,Y]=meshgrid(x,y);
load('LP_modes.mat')
OAM(:,:,1) = phi1(:,:,1);
OAM(:,:,2) = 1/sqrt(2)*(phi1(:,:,2)+1i*phi1(:,:,3));
OAM(:,:,3) = 1/sqrt(2)*(phi1(:,:,2)-1i*phi1(:,:,3));
OAM(:,:,4) = 1/sqrt(2)*(phi1(:,:,4)+1i*phi1(:,:,5));
OAM(:,:,5) = 1/sqrt(2)*(phi1(:,:,4)-1i*phi1(:,:,5));
OAM(:,:,6) = phi1(:,:,6);
phi = OAM(:,:,i);
I=sqrt(phi.*conj(phi));
ph=angle(phi);
v=abs(phi)/max(max(abs(phi)));%%LG field normalized to unity 
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
    Hol0=F.*sin(ph+2*pi*X*gx);
    Hol=Hol0-min(Hol0(:));
    SLM0=mod(Hol,1.17*pi);
    SLM=SLM0/max(SLM0(:))*150;
    

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