function [] = GetBodeFreq(Uin,tU,Yout,tY,uFreq)
%% GetBodeFreq: 
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
U.Time = t_p;
U.Pos = panel2deg(data(:,2));  % pattern x-pos: subtract mean and convert to deg [deg]  
U.Pos = FitPanel(U.Pos,U.Time ,tU); % fit panel data
U.Time = tU;
Y.Time = t_v;
Y.Pos = hAngles;
uFreq = 2;
nFreq = length(uFreq);
fTol = 0.2;
[U.Freq , U.Mag , U.Phase] = FFT(U.Time,U.Pos);
[Y.Freq , Y.Mag , Y.Phase] = FFT(Y.Time,Y.Pos);

figure ; clf ; hold on
plot(U.Time,U.Pos)
plot(Y.Time,Y.Pos)

[U.MAG,U.PHASE] = Get_IO_Freq(U.Freq,U.Mag,U.Phase,uFreq,0.2,true);
[Y.MAG,Y.PHASE] = Get_IO_Freq(Y.Freq,Y.Mag,Y.Phase,uFreq,0.2,true);

% plot(Y.Freq,Y.Mag)

% [U.pks,U.locs] = findpeaks(U.Mag);
% plot(U.Freq(U.locs),U.pks,'*g')
% fIdx = find(U.pks==max(U.pks));
% plot(U.Freq(U.locs(fIdx)),U.pks(fIdx),'ro')
% 
% [Y.pks,Y.locs] = findpeaks(Y.Mag);
% plot(Y.Freq(Y.locs),Y.pks,'*g')
% fIdx = find(Y.pks==max(Y.pks));
% plot(Y.Freq(Y.locs(fIdx)),Y.pks(fIdx),'ro')
end
















