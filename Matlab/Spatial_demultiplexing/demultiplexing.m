n = 6

SLM = cell(6);
SLM{1} = Laguerre_Gaussian(0,0);
SLM{2} = Laguerre_Gaussian(0,1);
SLM{3} = Laguerre_Gaussian(0,-1);
SLM{4} = Laguerre_Gaussian(0,2);
SLM{5} = Laguerre_Gaussian(0,-2);
SLM{6} = Laguerre_Gaussian(1,0);
SLM{7} = Laguerre_Gaussian(1,1);
SLM{8} = Laguerre_Gaussian(1,-1);
SLM{9} = Laguerre_Gaussian(1,2);
SLM{10} = Laguerre_Gaussian(1,-2);

modes = cell(6);
modes{1} = build_LG(0,0);
modes{2} = build_LG(0,1);
modes{3} = build_LG(0,-1);
modes{4} = build_LG(0,2);
modes{5} = build_LG(0,-2);
modes{6} = build_LG(1,0);
modes{7} = build_LG(1,1);
modes{8} = build_LG(1,-1);
modes{9} = build_LG(1,2);
modes{10} = build_LG(1,-2);

figure(1)
for  i = 1:n
    subplot(2,n,i)
    imagesc(abs(modes{i}.^2))
end

for  i = 1:n
    subplot(2,n,i+n)
    imagesc(SLM{i})
end

%%

load('LP_modes.mat')
OAM(:,:,1) = phi1(:,:,1);
OAM(:,:,2) = 1/sqrt(2)*(phi1(:,:,2)+1i*phi1(:,:,3));
OAM(:,:,3) = 1/sqrt(2)*(phi1(:,:,2)-1i*phi1(:,:,3));
OAM(:,:,4) = 1/sqrt(2)*(phi1(:,:,4)+1i*phi1(:,:,5));
OAM(:,:,5) = 1/sqrt(2)*(phi1(:,:,4)-1i*phi1(:,:,5));
OAM(:,:,6) = phi1(:,:,6);

SLM = cell(1,6);
SLM{1} = Linearly_polarized(1);
SLM{2} = Linearly_polarized(2);
SLM{3} = Linearly_polarized(3);
SLM{4} = Linearly_polarized(4);
SLM{5} = Linearly_polarized(5);
SLM{6} = Linearly_polarized(6);

modes = cell(1,6);
modes{1} = OAM(:,:,1);
modes{2} = OAM(:,:,2);
modes{3} = OAM(:,:,3);
modes{4} = OAM(:,:,4);
modes{5} = OAM(:,:,5);
modes{6} = OAM(:,:,6);

figure(1)
for  i = 1:6
    subplot(2,6,i)
    imagesc(abs(modes{i}.^2))
end

for  i = 1:6
    subplot(2,6,i+6)
    imagesc(SLM{i})
end

%%
error = linspace(0,2*pi,20);
for k = 1:20
i = 1;j = 1;
far = fftshift(fft2(modes{i}.*SLM{j}));
s = size(far);
coupling = zeros(n,n);
x = [s(2)/2+1-50-21:s(2)/2+1-50+21];
y = [s(1)/2+1-21:s(1)/2+1+21];
% x = [s(2)/2+1-21:s(2)/2+1+21];
% y = [s(1)/2+1-21:s(1)/2+1+21];
filter1 = fftshift(fft2(modes{1}.*exp(1i*SLM{1})));
filter = filter1(y,x);
filter = filter/sqrt(sum(sum(abs(filter.^2))));

for i = 1:n
    for j = 1:n
far = fftshift(fft2(modes{i}.*exp(1i*SLM{j}).*exp(1i*error(k)*rand(s))));
s = size(far);
% far((s(1)/2-25):(s(1)/2+25),(s(2)/2-25):(s(2)/2+25))=zeros(51);
I = abs(far.^2);
I = I/sum(sum(I));

% for k =-1:1
%     for l = -1:1
%        coupling(i,j) = coupling(i,j)+ I(s(1)/2+1+k,s(2)/2-50+1+l);
%     end 
% end    

dif_1 = far(y,x)/sqrt(sum(sum(abs(far.^2))));
coupling(i,j) = abs(sum(sum(filter.*dif_1)))^2;
% coupling(i,j) = I(s(1)/2+1,s(2)/2-50+1);
% figure(5)
% subplot(n,n,i+n*(j-1))
% 
% imagesc((I(y,x)))
% % xlim([s(2)/2+1-50-4,s(2)/2+1-50+4])
% % ylim([s(1)/2+1-4,s(1)/2+1+4])



    end
end
% 
% figure(7)
% imagesc((coupling))
% ylabel('#mode')
% xlabel('#mask')
% set(gca,'YDir','normal')
% colorbar

c(k) = trace(coupling)/6
e(k) = (sum(sum(coupling)) - trace(coupling))/30
end
figure(8)
subplot(1,2,1)
plot(error/pi,c)
ylabel('Avarege coupling')
xlabel('Fliker peak (pi rad)')
subplot(1,2,2)
plot(error/pi,e)
ylabel('Avarege error')
xlabel('Fliker peak (pi rad)')