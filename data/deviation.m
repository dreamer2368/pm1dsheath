clear all
close all
clc

spec = importdata('record.out');
Nt = spec(1); N = spec(2); Ng = spec(3); L = spec(4); dt = spec(5);
Ni = 202;

%B0
fileID = fopen('xp0.bin');
xp0 = fread(fileID,N*Nt,'double');
xp0 = reshape(xp0,[N,Nt]);

fileID = fopen('vp0.bin');
vp0 = fread(fileID,N*Nt,'double');
vp0 = reshape(vp0,[N,Nt]);

fileID = fopen('E0.bin');
E0 = fread(fileID,Ng*Nt,'double');
E0 = reshape(E0,[Ng,Nt]);

%B0 + dB0
fileID = fopen('xp.bin');
xp = fread(fileID,N*Nt,'double');
xp = reshape(xp,[N,Nt]);

fileID = fopen('vp.bin');
vp = fread(fileID,N*Nt,'double');
vp = reshape(vp,[N,Nt]);

fileID = fopen('E.bin');
E = fread(fileID,Ng*Nt,'double');
E = reshape(E,[Ng,Nt]);

%%
close all

xg = L/Ng*( (0:Ng-1) + 0.5 );
dE = abs(E0-E)/mean(mean(abs(E)));

% %video clip
% writerObj = VideoWriter('dE.avi');
% writerObj.FrameRate = 25;
% open(writerObj);

for i=1:Nt
%     figure(1)
% %     plot(xg,E0(:,i),'-k',xg,E(:,i),'-r');
%     plot(xg,abs(E0(:,i)-E(:,i)));
%     axis([0 L -5e-5 5e-5]);
    
%     figure(2)
%     semilogy(xg,dE(:,i),'.k','markersize',10);
%     axis([0 L 1e-16 5e-3]);
%     title('$\delta E/\overline{E}$ perturbation growth','Interpreter','latex');
%     xlabel('$x$','Interpreter','latex');
%     ylabel('$\delta E/\overline{E}$','Interpreter','Latex');
%     set(gca,'fontsize',25);
%     
% %     %videoclip
% %     frame = getframe(gcf);
% %     writeVideo(writerObj,frame);
%     pause(.001);
end

% % videoclip close
% close(writerObj);

%%
close all

dxp = min( abs(xp-xp0), L-abs(xp-xp0) )/L;
dvp = abs(vp-vp0)/sqrt(mean(mean(vp.^2)));

% %video clip
% writerObj = VideoWriter('dxpdvp.avi');
% writerObj.FrameRate = 25;
% open(writerObj);

for i=1:Nt
%     figure(3)
%     plot(xp0(:,i),vp0(:,i),'.k',xp(:,i),vp(:,i),'.r');
%     axis([0 L -.5 .5]);
    
%     figure(4)
%     loglog(dxp(:,i),dvp(:,i),'.k');
%     axis([1e-18 1e-2 1e-17 5e-2]);
%     title('$\delta x_p/L$, $\delta v_p/\overline{v}_p$ perturbation growth','Interpreter','Latex');
%     xlabel('$\delta x_p/L$','Interpreter','Latex');
%     ylabel('$\delta v_p/\overline{v}_p$','Interpreter','Latex');
%     set(gca,'fontsize',25);
%     xtick = [1e-15 1e-9 1e-3];
%     xticklabel = {'10^-^1^5','10^-^9','10^-^3'};
%     set(gca,'xtick',xtick,'xticklabel',xticklabel);
%     
% %     %videoclip
% %     frame = getframe(gcf);
% %     writeVideo(writerObj,frame);
%     pause(.01);
end

% % videoclip close
% close(writerObj);

%%
close all

T = dt*(1:Nt);
figure(5)
for i=5353
    semilogy(T,dxp(i,:),'--k',T,dvp(i,:),'-.r');
    hold on
end
xlabel('time');
ylabel('$\delta x_p$, $\delta v_p$','Interpreter','Latex');
set(gca,'fontsize',25);

%%
close all

Te = Nt*dt./log( dxp(:,Nt)./dxp(:,Ni) );
figure(6)
histogram(Te);
title('e-folding time of $\delta x_p$','Interpreter','Latex');
xlabel('$t_e$','Interpreter','Latex');
ylabel('number of particles');
set(gca,'fontsize',25);

Te = Nt*dt./log( dvp(:,Nt)./dvp(:,Ni) );
figure(7)
histogram(Te);
title('e-folding time of $\delta v_p$','Interpreter','Latex');
xlabel('$t_e$','Interpreter','Latex');
ylabel('number of particles');
set(gca,'fontsize',25);

Te = Nt*dt./log( dE(:,Nt)./dE(:,Ni) );
figure(8)
histogram(Te);
title('e-folding time of $\delta E$','Interpreter','Latex');
xlabel('$t_e$','Interpreter','Latex');
ylabel('number of particles');
set(gca,'fontsize',25);

%%
close all

T = dt*(1:Nt);
Momentum = squeeze( sum( vp, 1) );

figure(9)
plot(T,Momentum);