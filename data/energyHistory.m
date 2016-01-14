clear all
close all
clc

EnergyHistory = importdata('energy.out',' ',1);
T = EnergyHistory.data(:,1);
K = EnergyHistory.data(:,2);
P = EnergyHistory.data(:,3);
Total = K+P;

figure(1)
plot(T,Total,'-k',T,K,'-r',T,P,'-b');
title('Energy History','fontsize',20);
xlabel('Time','fontsize',20);
ylabel('Energy','fontsize',20);
legend('Total','Kinetic','Potential');