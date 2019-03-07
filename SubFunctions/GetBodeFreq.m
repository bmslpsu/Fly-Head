function [] = GetBodeFreq(U,tU,Y,tY,uFreq)
%% MakeData_Sine_HeadFree: Reads in all raw trials, transforms data, and saves in organized structure for use with figure functions
%   INPUTS:
%       U    	: input data
%       Y   	: output data
%   OUTPUTS:
%       PAT     : pattern structure
%       WING   	: wings structure
%---------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE INPUT %
%---------------------------------------------------------------------------------------------------------------------------------
close all
tU = (0:(1/200):10)';
U = panel2deg(data(:,2));  % pattern x-pos: subtract mean and convert to deg [deg]  
U = FitPanel(U,t_p,tU); % fit panel data
tY = t_v;
Y = hAngles;
[UFreq , UMag , UPhase] = FFT(tU,U);
[YFreq , YMag , YPhase] = FFT(tY,Y);

figure ; clf ; hold on
plot(tU,U)
plot(tY,Y)

figure ; clf ; hold on
plot(UFreq,UMag)
plot(YFreq,YMag)

[pks,locs] = findpeaks(UMag);
plot(UFreq(locs),pks,'*g')
fIdx = find(pks==max(pks));
plot(UFreq(locs(fIdx)),pks(fIdx),'ro')








end