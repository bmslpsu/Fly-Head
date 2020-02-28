function [acor,timeLag,maxCC,timeDiff] = CrossCorr(x,y,Fs)
%% CrossCorr: computes the cross-correlation between two signals and returns the correlations, lags, maximum corss-correlation,
% and time difference
%   INPUTS:
%       x       :   x-data
%       y       :   y-data
%       Fs   	:   sampling time [Hz]
%   OUTPUTS:
%       acor        :   cross-correlation
%       lag         :   time-lags
%       maxCC     	:   maximum cross-correlation
%       timeDiff 	:   time-difference
%

[acor,lag] = xcorr(x,y);
timeLag = lag/Fs;
[~,idx] = max(abs(acor));
maxCC = acor(idx);
timeDiff = timeLag(idx);
end