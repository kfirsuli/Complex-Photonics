% This method propagates 2 particle amplitudes in 1D random walk

n = 128; %system size
v = 1; %record video
%initiate a,b
phi_1 = zeros(n);
phi_1(62,63) = 1;
phi_1(63,62) = 1;
phi_1 = reshape(phi_1,n^2,1);

phi_2 = zeros(n);
phi_2(60,60) = 1;
phi_2(61,61) = +1;
phi_2(62,62) = 1;
phi_2(63,63) = +1;
phi_2(64,64) = 1;
phi_2(65,65) = +1;
phi_2(66,66) = 1;
phi_2(67,67) = +1;
phi_2 = reshape(phi_2,n^2,1);

phi_3 = zeros(n);
phi_3(62,62) = 1;
phi_3(64,64) = 1;
phi_3 = reshape(phi_3,n^2,1);

% propagator
I = speye(n,n);
E = sparse(2:n,1:n-1,1,n,n);
D = 1i*E+1i*E'+100*I;
A = kron(D,I)*kron(I,D);
% imagesc(abs(A))

if v == 1
vidfile = VideoWriter('testmovie.mp4','MPEG-4');
open(vidfile);
end
for i = 1:2000
    phi_1 = A*phi_1;
    phi_1 = phi_1/sqrt(sum(sum(abs(phi_1).^2)));

    
    phi_2 = A*phi_2;
    phi_2 = phi_2/sqrt(sum(sum(abs(phi_2).^2)));
    
    phi_3 = A*phi_3;
    phi_3 = phi_3/sqrt(sum(sum(abs(phi_3).^2)));    
    
    if mod(i,20) ==1 
    phi_1 = reshape(phi_1,n,n);
    phi_2 = reshape(phi_2,n,n);
    phi_3 = reshape(phi_3,n,n);
    
    subplot(1,3,1)
    imagesc(abs(phi_1).^2);
    set(gca,'YDir','normal')
    ylabel('Signal mode')
    xlabel('Idler mode')
    title({'2 distinguishable photons in 1D array', '\psi_{init} = |62,63>+|63,62>'})
    grid on
    
    subplot(1,3,2)
    imagesc(abs(phi_2).^2);
    set(gca,'YDir','normal')
    ylabel('Signal mode')
    xlabel('Idler mode')
    title({'2 distinguishable photons in 1D array', '\psi_{init} = |60,60>+|61,61>+...+|67,67>'})
    grid on
    
    subplot(1,3,3)
    imagesc(abs(phi_3).^2);
    set(gca,'YDir','normal')
    ylabel('Signal mode')
    xlabel('Idler mode')
    title({'2 distinguishable photons in 1D array', '\psi_{init} = |62>|64>+|64>|62>'})
    grid on
    
    if v == 1
    drawnow
    F(i) = getframe(gcf); 
    writeVideo(vidfile,F(i));
    end
    phi_1 = reshape(phi_1,n^2,1);
    phi_2 = reshape(phi_2,n^2,1);
    phi_3 = reshape(phi_3,n^2,1);
    pause(0.001)
    end
end
if v == 1
close(vidfile)
end
