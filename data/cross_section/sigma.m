close all
clear all
clc

N = 10000;
fileID = fopen('sig1.bin');
sig1 = fread(fileID,N,'double');
fileID = fopen('sig2.bin');
sig2 = fread(fileID,N,'double');
fileID = fopen('sig3.bin');
sig3 = fread(fileID,N,'double');
fileID = fopen('sig4.bin');
sig4 = fread(fileID,N,'double');
fileID = fopen('sig5.bin');
sig5 = fread(fileID,N,'double');
fileID = fopen('energy.bin');
energy = fread(fileID,N,'double');

%%
close all

figure(1)
loglog(energy,sig1*10^(20),'linewidth',2);
hold on
loglog(energy,sig2*10^(20),'linewidth',2);
loglog(energy,sig3*10^(20),'linewidth',2);
loglog(energy,sig4*10^(20),'linewidth',2);
loglog(energy,sig5*10^(20),'linewidth',2);
axis([0.01 1000 1e-3 1e2]);
xlabel('Incident particle energy(e^- or Ar^+)');
ylabel('$\sigma$($\times10^{-16}cm^2$)','interpreter','latex');
title('Cross section for collisions');
legend('e + Ar (elastic)','e + Ar (excite)','e + Ar (ionize)','Ar + Ar^+ (exchange)','Ar + Ar^+ (elastic)');
set(gca,'fontsize',25);