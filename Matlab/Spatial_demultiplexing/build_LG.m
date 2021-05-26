% clear all; close all; clc
function LG = build_LG(p,l)
 H=1024;
 V=1024;
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
La=Laguerre(p,abs(l),2*rho.^2/(w^2));
LG=w0/w*sqrt(2*factorial(p)/(pi*factorial(abs(l)+p)))*(2*rho.^2/w^2).^(abs(l)/2).*La.*exp(1i*(2*p+l+1)*atan(z/zr)).*exp(-rho.^2/w^2).*exp(-1i*k*rho.^2/R).*exp(-1i*l*phi);    
LG=LG/sqrt(sum(sum(abs(LG).^2)));
% imagesc(abs(LG.^2))