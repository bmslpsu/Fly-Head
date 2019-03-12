function [acor,timeLag,timeDiff] = CrossCorr(x,y,Fs)
%% CrossCorr: computes the corss-correlation between to signals and returns the the correlations, lags, and time difference
%   INPUTS:
%       x       :   x-data
%       y       :   y-data
%       Fs   	:   sampling time [Hz]
%   OUTPUTS:
%       acor    :   cross-correlation
%       lag   	:   time-lags
%       tdiff 	:   time-difference
%---------------------------------------------------------------------------------------------------------------------------------
[acor,lag] = xcorr(x,y);
timeLag = lag/Fs;
[~,idx] = max(abs(acor));
timeDiff = timeLag(idx);
end