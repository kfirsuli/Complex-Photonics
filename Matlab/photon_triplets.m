hold on
lambda = 532e-9;
b = 4e-4; % parabolic parameter
D = 2e-6; % diameter (https://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=358)
k = 2*pi/lambda;
n_0 = 1.46+0.06; % central refractive index 
N = 128;
n_count_modes = 15;
alpha = k*n_0/b;
sig = 1; % R(+1) or L(-1) circular pol.

img_size = 1*D;
range = linspace(-img_size/2, img_size/2, N);
[x, y] = meshgrid(range);
refractive_idx = n_0*(1 - range.^2/b^2);
[theta, rho] = cart2pol(x,y);

lmax = 20;
mmax = 15;
[l_grid, m_grid] = meshgrid(-lmax:lmax, 0:mmax);
lm_map = zeros(3, numel(l_grid));
[lm_map(1, :), I] = sort( abs(l_grid(:)) + 2*m_grid(:), 'ascend' ); % group by asending |l| + 2*m
lm_map(2, :) = l_grid(I);
lm_map(3, :) = m_grid(I);
lm_map = lm_map(:, 1:n_count_modes);

%% calculate the scalar field of each mode
n_modes = size(lm_map, 2);

% for ii = 1:n_modes
%     l = lm_map(2, ii); 
%     m = lm_map(3, ii); 
%     
%     E(:,:,ii) = GRIN_MMF_LG_fields( alpha, l, m, rho, theta );
% end
% betas = GRIN_MMF_LG_betas( lambda, n_0, b, alpha, lm_map, sig, inf ); % propagation constants for straight MMF
% E19 = GRIN_MMF_LG_fields( alpha, 0, 19, rho, theta );
% % betas19 = GRIN_MMF_LG_betas( lambda, n_0, b, alpha, 741, sig, inf ); % propagation constants for straight MMF
% E19 = E_normalize(E19);
% E2 = GRIN_MMF_LG_fields( alpha, 0, 1, rho, theta );
% % betas19 = GRIN_MMF_LG_betas( lambda, n_0, b, alpha, 741, sig, inf ); % propagation constants for straight MMF
% E2 = E_normalize(E2);
% E1 = GRIN_MMF_LG_fields( alpha/3, 0, 0, rho, theta );
% % betas1 = GRIN_MMF_LG_betas( lambda, n_0*0.99, b, alpha*0.99/3, 1, sig, inf ); % propagation constants for straight MMF
% E1 = E_normalize(E1);

b1 = 1/pi/2*lambda*GRIN_MMF_LG_betas(lambda, n_0, b, alpha, lm_map, sig, 0);
figure(2)
plot(b1)


% figure(3)
% subplot(1,3,1)
% title('LP019 @ 532nm')
% imagesc(range, range,abs(E19.^2))
% subplot(1,3,2)
% title('LP02 @ 532nm')
% imagesc(range, range,abs(E2.^2))
% subplot(1,3,3)
% title('LP01 @ 532nm')
% imagesc(range, range,abs(E1.^2))
% 
% sum(sum(E1*conj(E1)))
% sum(sum(E19*conj(E19)))
% sum(sum(E2*conj(E2)))
% sum(sum(E19.*conj(E1).*conj(E1).*conj(E1)))
% sum(sum(E2.*conj(E1).*conj(E1).*conj(E1)))


%% some self-defined functions
function E = GRIN_MMF_LG_fields( alpha, l, m, r, theta )
    ar_2 = alpha*(r.^2);

    c1 = sqrt( (alpha/(2*pi))*( (2*factorial(m))/factorial(m+abs(l)) ) );
    c2 = exp( -ar_2/2 );
    c3 = ar_2.^(abs(l)/2);
    c4 = laguerreL(m, abs(l), ar_2); % Generalized Laguerre polynomials (https://en.wikipedia.org/wiki/Laguerre_polynomials)
    c5 = exp(1i*l*theta);
    E = c1*c2.*c3.*c4.*c5; % scalar field
    E = E_normalize(E);
end

function [betas, betas_prime] = GRIN_MMF_LG_betas( lambda, n_0, b, alpha, lm_map, sig, bend_r )
% lm_map is a 3 by n_modes matrix, with each column as the characteristics of one particular mode
%                                      1st row as the angular momentum, 2m + |l|
%                                      2nd row as the orbital angular momentum
%                                      3rd row as the # of radial nodes
    k = 2*pi/lambda;
    l_s = lm_map(2,:);
    m_s = lm_map(3,:);
    
    betas = sqrt((k*n_0)^2 - 2*alpha*(abs(l_s) + 2*m_s + 1)); % 1 by n_modes matrix, propagation constants in straight MMF
    curv = 1/bend_r; % curvature
    SO_corrections = -(l_s*sig + 1)/(2*k*n_0*b^2); % 1 by n_modes matrix, spin-orbital interactions 
    betas_prime = betas + (curv^2)*( (k*n_0*b^2)/2 - 9*b*(abs(l_s) + 2*m_s + 1)/4 ) + SO_corrections; % 1 by n_modes matrix, propagation constants in bent MMF
end

function E = E_normalize(E)
    E = E/sqrt(sum(sum(E*conj(E))));
end

% %%
% c = 2.997e+8;
% 
% c_1 = 1.132e-6;
% c_2 = -2.059e-6; 
% c_4 = 6.734e-6; 
% c_5 = -2.239e-7;
% x = linspace(0.4,2.4,99);
% n_UHNA1 = c_1*x.^-4+c_2*x.^-2+0*x.^0+c_4*x.^2+c_5*x.^4;
% D = (1/c)*(-20*c_1*x.^-5-6*c_2*x.^-3-2*c_4*x.^1-12*c_5*x.^3);
% figure(10)
% subplot(1,2,1)
% hold on
% plot(x*1000,n_UHNA1,'Color',[0 0.4470 0.7410])
% plot(3*x(1:33)*1000,n_UHNA1(1:33),'--','Color',[0 0.4470 0.7410])
% 
% subplot(1,2,2)
% hold on
% plot(x*1000,D*1e+15)
% 
% % c_1 = -5.855e-6;
% % c_2 = 17.273e-6; 
% % c_4 = 5.029e-6; 
% % c_5 = -1.186e-7;
% % x = linspace(0.4,2.4,100);
% % n_UHNA3 = c_1*x.^-4+c_2*x.^-2+0*x.^0+c_4*x.^2+c_5*x.^4;
% % D = (1/c)*(-20*c_1*x.^-5-6*c_2*x.^-3-2*c_4*x.^1-12*c_5*x.^3);
% % figure(10)
% % subplot(1,2,1)
% % plot(x*1000,n_UHNA3)
% % subplot(1,2,2)
% % plot(x*1000,D*1e+15)
% 
% c_1 = 1.982e-6;
% c_2 = -5.010e-6; 
% c_4 = 10.433e-6; 
% c_5 = -2.798e-7;
% x = linspace(0.4,2.4,100);
% n_UHNA4 = c_1*x.^-4+c_2*x.^-2+0*x.^0+c_4*x.^2+c_5*x.^4;
% D = (1/c)*(-20*c_1*x.^-5-6*c_2*x.^-3-2*c_4*x.^1-12*c_5*x.^3);
% figure(10)
% subplot(1,2,1)
% plot(x*1000,n_UHNA4,'Color',[0.8500 0.3250 0.0980])
% plot(3*x(1:33)*1000,n_UHNA4(1:33),'--','Color',[0.8500 0.3250 0.0980])
% 
% subplot(1,2,2)
% plot(x*1000,D*1e+15)
% 
% 
% a1=0.6961663;
% a2=0.4079426;
% a3=0.8974794;
% b1= 0.0684043;
% b2=0.1162414;
% b3=9.896161;
% 
% subplot(1,2,2)
% ylabel('D (ps/nm/km)')
% xlabel('Wavelength (nm)')
% subplot(1,2,1)
% ylabel('\Deltan')
% xlabel('Wavelength (nm)')