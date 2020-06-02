clear;close all; clc

Fs = 500;
T = 5;
t = (0:(1/Fs):T)';
f = 1;
phi = deg2rad(-30);
td = (phi/(2*pi)) * (1/f);
x = 1*sin(2*pi*f*t + 0);
y = 1*sin(2*pi*f*t + phi);

[acor,timeLag,maxCC,timeDiff] = CrossCorr(x, y, Fs,'normalized');

disp(timeDiff)

figure (1) ; cla ; hold on
plot(t,x)
plot(t,y)