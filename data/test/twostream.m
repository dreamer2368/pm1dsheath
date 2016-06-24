clear all
close all
clc

spec = importdata('record');
N = spec(1); Ng = spec(2); Nt = spec(3); L = spec(4);

fileID = fopen('xp.bin');
xp = fread(fileID,N*Nt,'double');
xp = reshape(xp,[N,Nt]);

fileID = fopen('vp.bin');
vp = fread(fileID,N*Nt,'double');
vp = reshape(vp,[N,Nt]);

fileID = fopen('E.bin');
E = fread(fileID,Ng*Nt,'double');
E = reshape(E,[Ng,Nt]);

fileID = fopen('PE.bin');
PE = fread(fileID,Nt,'double');

fileID = fopen('KE.bin');
KE = fread(fileID,Nt,'double');

fileID = fopen('rho.bin');
rho = fread(fileID,Ng,'double');

%%
close all

for i=1:Nt
    plot(xp(:,i),vp(:,i),'.k');
    pause(.1);
end