function [newFv,newFREQ,N,bins] = DFT_bin(Fv,FREQ,bin_edges)
%% DFT_bin: averages DFT into bins
%
%   INPUT:
%       Fv              :   frequency vector
%       FREQ         	:   data to bin (complex or real)
%       bin_edges       :   frequency bin edges
%
%   OUTPUT:
%       newFv           :   new frequency bins
%       newFREQ       	:   new binned average of frequency domain data
%       N               :   bin counts
%       bins            :   bin locations
%

% freq_bin_edges = ((0:0.2:100)-0.001)';

[N, ~, bins] = histcounts(Fv,bin_edges);
N = N';

bins(end) = bins(end-1);

newFv   = splitapply(@mean, Fv,   bins);
newFREQ = splitapply(@mean, FREQ, bins);

end