clear all
close all
clc

B = linspace(0.9,1.1,1001);
Jdata = importdata('Jdata.out');

figure(2)
for i = 1:size(Jdata,2)
    h = plot(B,Jdata(:,i)-mean(Jdata(:,i)));
    if( i==size(Jdata,2)-1 )
        h.Color = 'b';
    elseif( i==size(Jdata,2) )
        h.Color = 'r';
    else
        h.Color = 'k';
    end
    hold on
end
legend('1/2\pi T_p','0.5T_p','T_p','2T_p','4T_p','6T_p','8T_p','10T_p');
xlabel('B','fontsize',20);
ylabel('deviation','fontsize',20);
title('Noise generation by chaos in time','fontsize',20);
