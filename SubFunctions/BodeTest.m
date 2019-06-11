clear;close all;clc

T = 20;
Fs = 200;
t = (0:(1/Fs):T)';

U.A = 2;
U.f = 1;
U.phi = 0;
U.x = U.A*sin(2*pi*U.f*t + U.phi);
[U.Fv, U.Mag, U.Phs, U.Freq] = FFT(t,U.x);
[U.MAG,U.PHASE,U.fIdx] = Get_IO_Freq(U.Fv, U.Mag, U.Phs, U.f,0.2,true);
U.FREQ = U.Freq(U.fIdx);
U.R = real(U.FREQ);
U.I = imag(U.FREQ);

Y.A = 6;
Y.f = 1;
Y.phi = -pi/4;
Y.x = Y.A*sin(2*pi*Y.f*t + Y.phi);
[Y.Fv, Y.Mag, Y.Phs, Y.Freq] = FFT(t,Y.x);
[Y.MAG,Y.PHASE,Y.fIdx] = Get_IO_Freq(Y.Fv, Y.Mag, Y.Phs, Y.f);
Y.FREQ = Y.Freq(Y.fIdx);
Y.R = real(Y.FREQ);
Y.I = imag(Y.FREQ);

Gain = Y.MAG./U.MAG
Phase = rad2deg(Y.PHASE - U.PHASE)

CmplxGain = Y.FREQ./U.FREQ;
cGain = abs(CmplxGain)
cPhase = rad2deg(angle(CmplxGain))
Fly
% rr(:,kk) = real(freq(:,kk));
% img(:,kk) = imag(freq(:,kk));
% [MAG(:,kk),PHASE(:,kk),fIdx(:,kk)] = Get_IO_Freq(Fv(:,kk),Mag(:,kk),Phs(:,kk),ff(kk));
% FREQ = freq(fIdx(:,kk),:);
% R(:,kk) = real(FREQ(:,kk));
% I(:,kk) = imag(FREQ(:,kk));