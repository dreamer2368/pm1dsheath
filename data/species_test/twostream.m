clear all
close all
clc

spec = importdata('record');
N = 5000; Ng = spec(2); Nt = spec(3); L = spec(4);

fileID = fopen('xp_1.bin');
xp1 = fread(fileID,N*Nt,'double');
xp1 = reshape(xp1,[N,Nt]);

fileID = fopen('vp_1.bin');
vp1 = fread(fileID,N*Nt,'double');
vp1 = reshape(vp1,[N,Nt]);

fileID = fopen('xp_2.bin');
xp2 = fread(fileID,N*Nt,'double');
xp2 = reshape(xp2,[N,Nt]);

fileID = fopen('vp_2.bin');
vp2 = fread(fileID,N*Nt,'double');
vp2 = reshape(vp2,[N,Nt]);

fileID = fopen('E.bin');
E = fread(fileID,Ng*Nt,'double');
E = reshape(E,[Ng,Nt]);

fileID = fopen('PE.bin');
PE = fread(fileID,Nt,'double');

fileID = fopen('KE_1.bin');
KE1 = fread(fileID,Nt,'double');

fileID = fopen('rho.bin');
rho = fread(fileID,Ng,'double');

%%
close all

for i=1:Nt
    plot(xp1(:,i),vp1(:,i),'.k',xp2(:,i),vp2(:,i),'.r');
    pause(.1);
end