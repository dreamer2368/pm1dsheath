clear all
close all
clc

spec = importdata('record');
N = spec(1); Ng = spec(2); Nt = spec(3); L = spec(4); mod = spec(5);
Nt = floor(Nt/mod);

fileID = fopen('Np.bin');
Np = fread(fileID,N*Nt,'int32');
Np = reshape(Np, [N,Nt]);

fileID = fopen('E.bin');
E = fread(fileID,Ng*Nt,'double');
E = reshape(E,[Ng,Nt]);

fileID = fopen('PE.bin');
PE = fread(fileID,Nt,'double');

fileID = fopen('phi.bin');
phi = fread(fileID,Ng*Nt,'double');
phi = reshape(phi,[Ng,Nt]);

fileID = fopen('rho.bin');
rho = fread(fileID,Ng*Nt,'double');
rho = reshape(rho,[Ng,Nt]);

%%
close all
clc
npsum1 = 1; npsum2 = 1;
dx = L/(Ng-1);
xg = dx*(0:Ng-1);

%video clip
writerObj = VideoWriter('field.avi');
writerObj.FrameRate = 20;
open(writerObj);

for i=1:Nt
    figure(2)
    plot(xg,rho(:,i),'-k',xg,E(:,i),'-b',xg,phi(:,i),'-r');
%     axis([0 L -.5e-5 2.5e-5]);
    axis([0 L -3e-2 3e-2]);
    title('Field');
    xlabel('$x$(m)','interpreter','latex');
    ylabel('$\rho$, $E$, $\phi$','interpreter','latex');
    set(gca,'fontsize',25);
    
    
    %videoclip
    frame = getframe(gcf);
    writeVideo(writerObj,frame);
%     pause();
end

% videoclip close
close(writerObj);

%%
close all
dx = L/(Ng-1);
xg = dx*(0:Ng-1);

i=floor(1);
    fileID = fopen(strcat('xp/',num2str(i),'_1.bin'));
    xp_e = fread(fileID,Np(1,i+1),'double');
    fileID = fopen(strcat('vp/',num2str(i),'_1.bin'));
    vp_e = fread(fileID,Np(1,i+1),'double');
    
    f1=figure(1);
    plot(xp_e,vp_e,'.k');
%     axis([0 L -2e-1 2e-1]);
    title('Electron distribution');
    xlabel('$x$(m)','interpreter','latex');
    ylabel('$v$(m/s)','interpreter','latex');
    set(gca,'fontsize',25);
    
    figure(2)
    histogram(vp_e);
    
    fclose('all');

%%
close all

dx = L/Ng;
NFFT = 2^nextpow2(Ng);
Ekhistory = [];
for i=1:Nt
    Ek = fft(E(:,i),NFFT)/Ng;
    Ekhistory = [Ekhistory 2*abs(Ek(1:NFFT/2+1))];
end
dt = 0.2; time = (1:Nt)*dt;
k = 2*pi*(1/dx)/2*linspace(0,1,NFFT/2+1);

figure(3)
semilogy(time, 2*(Ekhistory(2,:)),'-k',time(1:Nt/2),0.00001*exp(1/2/sqrt(2)*time(1:Nt/2)),'-r');
title('$\hat{E}(k=\frac{2pi}{L})$','Interpreter','Latex');
xlabel('time'); ylabel('$\hat{E}(k)$','Interpreter','Latex');
set(gca,'fontsize',25);