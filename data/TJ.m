close all
clear all
clc

T = importdata('T.out');
J0 = importdata('J0.out');
J1 = importdata('J1.out');
absDJDA = abs((J1-J0)/1e-15);
figure(1)
plot(T,J0,T,J1);
title('J in time','fontsize',20);
xlabel('time','fontsize',20);
ylabel('J','fontsize',20);
legend('J0','J1');
figure(2)
semilogy(T,absDJDA);
title('|\Delta J/\Delta A|','fontsize',20);
xlabel('time','fontsize',20);
ylabel('|\Delta J/\Delta A|','fontsize',20);