clear all
close all
clc

spec = importdata('record');
N = spec(1); Ng = spec(2); Nt = spec(3); L = spec(4); mod = spec(5);
Nt = floor(Nt/mod);

fileID = fopen('Np.bin');
Np = fread(fileID,N*Nt,'int32');
Np = reshape(Np, [N,Nt]);

% fileID = fopen('xp_1.bin');
% xp1 = fread(fileID,sum(Np(1,:)),'double');
% 
% fileID = fopen('vp_1.bin');
% vp1 = fread(fileID,sum(Np(1,:)),'double');
% 
% fileID = fopen('xp_2.bin');
% xp2 = fread(fileID,sum(Np(2,:)),'double');
% 
% fileID = fopen('vp_2.bin');
% vp2 = fread(fileID,sum(Np(2,:)),'double');

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
dx = L/Ng;
xg = dx*(1:Ng);

% %video clip
% writerObj = VideoWriter('phi.avi');
% writerObj.FrameRate = 20;
% open(writerObj);

for i=1:Nt
%     figure(1)
%     plot(xp1(npsum1:npsum1-1+Np(1,i)),vp1(npsum1:npsum1-1+Np(1,i)),'.k',xp2(npsum2:npsum2-1+Np(2,i)),vp2(npsum2:npsum2-1+Np(2,i)),'.r');
%     npsum1 = npsum1 + Np(1,i);
%     npsum2 = npsum2 + Np(2,i);
    
    figure(2)
    plot(xg,rho(:,i),'-k');
    axis([0 L -1e-5 2e-5]);
%     axis([0 L 0 4]);
    title('charge density');
    xlabel('$x$(m)','interpreter','latex');
    ylabel('$\rho$(C/m)','interpreter','latex');
    set(gca,'fontsize',25);
    
    
%     %videoclip
%     frame = getframe(gcf);
%     writeVideo(writerObj,frame);
    pause();
end

% % videoclip close
% close(writerObj);

%%
close all

figure(1)
plot(phi(200,:));

%%
close all
ve = 593097.26001485332; vi = 1553.1051663509243;

% %video clip
% writerObj1 = VideoWriter('electron.avi');
% writerObj1.FrameRate = 20;
% open(writerObj1);
% 
% %video clip
% writerObj2 = VideoWriter('ion.avi');
% writerObj2.FrameRate = 20;
% open(writerObj2);

for i=1:Nt
    fileID = fopen(strcat('xp/',num2str(i),'_1.bin'));
    xp_e = fread(fileID,Np(1,i),'double');
    fileID = fopen(strcat('vp/',num2str(i),'_1.bin'));
    vp_e = fread(fileID,Np(1,i),'double');
    
    fileID = fopen(strcat('xp/',num2str(i),'_2.bin'));
    xp_i = fread(fileID,Np(2,i),'double');
    fileID = fopen(strcat('vp/',num2str(i),'_2.bin'));
    vp_i = fread(fileID,Np(2,i),'double');
    
    f1=figure(1);
    plot(xp_e,vp_e,'.k');
    axis([0 L -3*ve 3*ve]);
    title('Electron distribution');
    xlabel('$x$(m)','interpreter','latex');
    ylabel('$v$(m/s)','interpreter','latex');
    set(gca,'fontsize',25);
    
%     %videoclip
%     frame = getframe(f1);
%     writeVideo(writerObj1,frame);
    
    f2=figure(2);
    plot(xp_i,vp_i,'.r');
    axis([0 L -5*vi 5*vi]);
    title('Ion distribution');
    xlabel('$x$(m)','interpreter','latex');
    ylabel('$v$(m/s)','interpreter','latex');
    set(gca,'fontsize',25);
    
%     %videoclip
%     frame = getframe(f2);
%     writeVideo(writerObj2,frame);
    
    fclose('all');
    pause(.1);
end

% % videoclip close
% close(writerObj1);
% close(writerObj2);

%%
close all
fileID = fopen(strcat('vp/1_1.bin'));
vp_e = fread(fileID,Np(1,1),'double');

figure(1)
histogram(vp_e);