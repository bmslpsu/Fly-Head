function [PeakPatMag, PeakPatFreq, PeakHeadMag, PeakHeadFreq, PeakWingMag, PeakWingFreq] = SOSPeakFinder(npeaks, patFreq, patMag, headFreq, headMag, wingFreq, wingMag)
% SOSPeakFinder Input number of pattern frequencies and FFT data of Pattern,
% Head and Wings. Include Frequency then Magnitude.
% It finds the peak within +/- .3 of the peak pattern frequency


[PeakPatMag,PeakPatLoc] = findpeaks(patMag, 'NPeaks',npeaks, 'SortStr','descend');
PeakPatFreq = patFreq(PeakPatLoc);

PeakHeadMag = zeros(npeaks, 1);
PeakHeadFreq = zeros(npeaks, 1);
PeakWingMag = zeros(npeaks, 1);
PeakWingFreq = zeros(npeaks, 1);

for kk = 1:npeaks
    x = PeakPatFreq(kk);
    xlow = x-.3;
    xhigh = x+.3;
    loc1 = find(headFreq >= xlow  & headFreq<= xhigh);
    PeakHeadMag(kk,1) = max(headMag(loc1));
    ind1 = find(headMag == PeakHeadMag(1,kk));
    PeakHeadFreq(kk,1) = flyFreq(ind1);
    loc2 = find(wingFreq >= xlow  & wingFreq<= xhigh);
    PeakWingMag(kk,1) = max(wingdMag(loc2));
    ind2 = find(wingMag == PeakWingMag(1,kk));
    PeakWingFreq(kk,1) = wingFreq(ind2);
end


hold on
plot(patFreq,patMag,PeakPatFreq,PeakPatMag,'or')
plot(headFreq, headMag, PeakHeadFreq,PeakHeadMag, 'o')
plot(wingFreq, wingMag, PeakWingFreq,PeakWingMag, 'o')
hold off

end