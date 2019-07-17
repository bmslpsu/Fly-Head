function [EMD_ALL_filt] = HR_sim(Eye_ALL)

% TEMPORALLY FILTER SPACE-TIME EYE SIGNAL
t  = 0:0.001:0.060;
tp = 0.020; % 20 ms time to peak; from Dickson et al. 2008
sigma = 0.355; % width
% linear impulse response light energy -> neural signal
Vt = exp(-((log(t./tp).^2) ./ (2*sigma^2)));

% convolve signal for each ommatidium
[~, num_receptors] = size(Eye_ALL);
eye_filt = nan(size(Eye_ALL));

for jj = 1:num_receptors % for each eye
    eye_filt(:,jj) = conv(Eye_ALL(:,jj), Vt, 'same');
end

% COMPUTE EMD RESPONSE

% initializations for HR model
samp_rate = 1000;
tau = 0.035;
h = 1/samp_rate;
% discrete low-pass filter parameters
a = h / (tau+h);

% initial conditions
InMat = 5*(rand(1,num_receptors) - 0.5);
FiltMat = zeros(size(InMat)); FiltMat_1 = zeros(size(InMat));

EMD_ALL_filt = nan(length(eye_filt), num_receptors);

for kk = 1:length(eye_filt) % for each time point
    
    % compute HR motion
    InMat = eye_filt(kk, :);
    FiltMat = a*(InMat) + (1-a)*FiltMat_1; % discrete low-pass filter
    FiltMat_1 = FiltMat;
    HR_Motion = Va(t-tau) * Vb - Vb(t-tau) * Va;
%     HR_Motion = (FiltMat(1:end-1).*InMat(2:end) - FiltMat(2:end).*InMat(1:end-1)); % delay and correlate
    
    EMD_ALL_filt(kk,1:end-1) = HR_Motion; 
    
end

end