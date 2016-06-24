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
dx = L/(Ng-1);
xg = dx*(0:Ng-1);

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
    plot(xg,phi(:,i) - phi(floor(Ng/2),i),'-k');
%     axis([0 L -1e0 1]);
%     axis([0 L 0 4]);
    title('potential');
    xlabel('$x$(m)','interpreter','latex');
    ylabel('$\phi$(V)','interpreter','latex');
    set(gca,'fontsize',25);
    
    
%     %videoclip
%     frame = getframe(gcf);
%     writeVideo(writerObj,frame);
    pause(.0001);
end

% % videoclip close
% close(writerObj);

%%
close all

figure(1)
plot(phi(200,:));

%%
close all
clc
dx = L/(Ng-1);
xg = dx*(0:Ng-1);

EV_TO_K = 11604.52;
Te = 50.0*EV_TO_K;
tau = 100.0;
Ti = Te/tau;
me = 9.10938215E-31;
mu = 1836;
mi = mu*me;
K = 1.38065E-23;
vB = sqrt(K*(Te+3*Ti)/mi);

ve = sqrt(K*Te/me); vi = sqrt(K*Ti/mi);

%video clip
writerObj1 = VideoWriter('electron.avi');
writerObj1.FrameRate = 60;
open(writerObj1);

%video clip
writerObj2 = VideoWriter('ion.avi');
writerObj2.FrameRate = 60;
open(writerObj2);

%video clip
writerObj3 = VideoWriter('phi.avi');
writerObj3.FrameRate = 60;
open(writerObj3);

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
    axis([0 L -5*ve 5*ve]);
    title('Electron distribution');
    xlabel('$x$(m)','interpreter','latex');
    ylabel('$v$(m/s)','interpreter','latex');
    set(gca,'fontsize',25);
    
    f2=figure(2);
    plot(xp_i,vp_i,'.r',[0 L], [vB vB], '-b', [0 L], [-vB -vB],'-b');
    axis([0 L -35*vi 35*vi]);
    title('Ion distribution');
    xlabel('$x$(m)','interpreter','latex');
    ylabel('$v$(m/s)','interpreter','latex');
    set(gca,'fontsize',25);
    
    f3=figure(3);
    plot(xg,phi(:,i),'-k');
    axis([0 L -2e2 1]);
    title('potential');
    xlabel('$x$(m)','interpreter','latex');
    ylabel('$\phi$(V)','interpreter','latex');
    set(gca,'fontsize',25);
    
    %videoclip
    frame = getframe(f1);
    writeVideo(writerObj1,frame);
    
    %videoclip
    frame = getframe(f2);
    writeVideo(writerObj2,frame);

    %videoclip
    frame = getframe(f3);
    writeVideo(writerObj3,frame);

    fclose('all');
%     pause(.00000001);
end

% videoclip close
close(writerObj1);
close(writerObj2);
close(writerObj3);

%%
close all

figure(1)
plot(1:Nt,Np(1,:),'-k',1:Nt,Np(2,:),'-r','linewidth',2);
xlabel('Time step');
ylabel('Particles');
title('Number of ion/electron');
legend('Electron','Ion');
set(gca,'fontsize',25);
figure(2)
plot(1:Nt,abs(Np(1,:)-Np(2,:)),'-k');
figure(3)
plot(1:Nt,phi(20,:),'-k');
axis([0 Nt 0 25]);

%%
close all

figure(1)
sample = ( abs(xp_i-0.5*L)<0.001 );
histogram( xp_i );