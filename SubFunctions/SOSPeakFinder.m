function [VarPeakFreq, VarPeakMag, VarPeakPhase] = SOSPeakFinder(npeaks, PatPeakFreq, VarFreq, VarMag, VarPhase)
% SOSPeakFinder Input pattern frequencies and FFT data of a single parameter Frequency then Magnitude.
% It finds the peak within +/- .3 of the peak pattern frequency

VarPeakMag = zeros(npeaks, 1);
VarPeakFreq = zeros(npeaks, 1);
VarPeakPhase = zeros(npeaks,1);


for kk = 1:npeaks
    x = PatPeakFreq(kk);
    xlow = x-.1;
    xhigh = x+.1;
    loc1 = find(VarFreq >= xlow  & VarFreq<= xhigh);
    VarPeakMag(kk,1) = max(VarMag(loc1));
    ind1 = find(VarMag == VarPeakMag(kk,1));
    VarPeakFreq(kk,1) = VarFreq(ind1);
    VarPeakPhase(kk,1) = VarPhase(ind1);
    
end




end