clear all
close all
clc

fDA = importdata('fDA.out');
ek = importdata('ek_twostream.out');
% ek = importdata('ek.out');
% ek = reshape(ek,[length(ek)/length(fDA), length(fDA)]);
Time = importdata('Tf.out');

Fig = figure(1);
set(Fig,'Position',[100,100,650,550]);
loglog(fDA,.0001*fDA,'-r','LineWidth',5);
hold on
% loglog(fDA,ek0,'.-k');
loglog(fDA,ek(1,:),'+-','LineWidth',5);
loglog(fDA,ek(length(Time),:),'o-','LineWidth',5);
for i=2:length(Time)-1;
    loglog(fDA,ek(i,:),'--');
end
loglog(fDA(length(fDA)),ek(length(Time),length(fDA)),'Or','markersize',25,'LineWidth',5);
axis([1e-10 1e5 1e-15 1e5]);
title('TSC, $N=10000$, $v_0=0.2$','Interpreter','Latex');
xlabel('$\Delta$','Interpreter','Latex');
ylabel('error','Interpreter','Latex');
% h=legend('$\mathcal{O}(\Delta)$','$T_p/2\pi$','$10T_p$');
h=legend('$\mathcal{O}(\Delta B)$','$T_p/2\pi$','$20T_p$','$0.5T_p$','$T_p$','$2T_p$','$4T_p$','$6T_p$','$8T_p$','$10T_p$');
set(h,'fontsize',15,'Interpreter','Latex');
set(gca,'fontsize',35);