function [EMD_ALL_filt] = HR(Eye_ALL,tt)
% HR_sim: 
%   INPUTS:
%       Eye_ALL         :   spatially filtered intesnity gradients from omnatidia
%       tt              :   time vector
%   OUTPUTS:
%       EMD_ALL_filt 	:   EMD output
%---------------------------------------------------------------------------------------------------------------------------------
% TEMPORALLY FILTER SPACE-TIME EYE SIGNAL
tp = 0.020; % 20 ms time to peak; from Dickson et al. 2008
sigma = 0.355; % width

% Linear impulse response light energy -> neural signal
Vt = exp(-((log(tt./tp).^2) ./ (2*sigma^2)));

% convolve signal for each ommatidium
[~, num_receptors] = size(Eye_ALL);
eye_filt = nan(size(Eye_ALL));

for jj = 1:num_receptors % for each eye
    eye_filt(:,jj) = conv(Eye_ALL(:,jj), Vt, 'same');
end

% COMPUTE EMD RESPONSE
% Initializations for HR model
samp_rate = 1/mean(diff(tt));
tau = 0.035;
h = 1/samp_rate;
a = h / (tau+h);

% Initial conditions
EMD_ALL_filt    = nan(length(eye_filt), num_receptors);
InMat           = 5*(rand(1,num_receptors) - 0.5);
% FiltMat         = zeros(size(InMat)); 
FiltMat_1       = zeros(size(InMat));
for kk = 1:size(eye_filt,1) % for each time point
    % Compute HR motion
    InMat = eye_filt(kk, :);
    FiltMat = a*(InMat) + (1-a)*FiltMat_1; % discrete low-pass filter
    FiltMat_1 = FiltMat;
    HR_Motion = (FiltMat(1:end-1).*InMat(2:end) - FiltMat(2:end).*InMat(1:end-1)); % delay and correlate
    
    EMD_ALL_filt(kk,1:end-1) = HR_Motion; 
end

end