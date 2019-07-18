function [EMD_ALL, Pattern_ALL, Eye_ALL] = EMDsim_v5_7_5deg_Grnd_new(pattern, vel, showplot, export)
% This simulation reconstructs the fly's visual environment and calculates
% optic flow based on Reichart detector

% Make eye filters
[eye_filt, num_samp_pts, num_receptors] = make_eye_filters();

% Initializations for HR model
samp_rate = 1000;
tau = 0.035;
h = 1/samp_rate;

% Discrete low-pass filter parameters
a = h / (tau+h);

InMat = 5*(rand(1,num_receptors) - 0.5);
FiltMat_1 = zeros(size(InMat));

% Specify position
step = vel/3.75;
Pos = ones(step*96,1); % initialize
for jj = 1:96-1
    Pos(jj*step:(jj+1)*step) = jj;
end

% Preallocate array
EMD_ALL = nan(length(Pos), num_receptors);
Pattern_ALL = nan(length(Pos), 96); % 96-pixel column
% Eye_ALL = nan(length(Pos), num_receptors);
Eye_ALL = nan(length(Pos), num_receptors);

mm = 1;
for jj = 1:length(Pos)
    
	if showplot == 1
        H = figure(50); cla; clf;
        set(H, 'Renderer','OpenGL');
        set(H, 'Position',[100, 100, 16*40, 16*50]);        
	end
	
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATE OUTPUT FROM ELEMENTARY MOTION DETECTOR
    % (HASSENTSEIN-REICHARDT)
    
    current_frame = squeeze(pattern.Pats(:, :, Pos(jj)));
    % Upsample by factor of 10
    for k = 1:10
        Up_frame(k:10:num_samp_pts) = current_frame(1,1:96);
    end
    
    % Get filtered eye projection
    eye_sample = Up_frame*eye_filt;
    
    % Compute HR motion
    InMat = eye_sample;
    FiltMat = a*(InMat) + (1-a)*FiltMat_1; % discrete low-pass filter
    FiltMat_1 = FiltMat;
    HR_Motion = (FiltMat(1:end-1).*InMat(2:end) - FiltMat(2:end).*InMat(1:end-1)); % delay and correlate

    % Prepare data for space-time plot
    EMD_ALL(jj,1:end-1) = HR_Motion;
    Pattern_ALL(jj,:) = current_frame(1,1:96);
    Eye_ALL(jj, :) = eye_sample;
    
    mm = mm + 1;
end

end

function [eye_filt, num_samp_pts, num_receptors] = make_eye_filters()
% rough approximation. follows from caption of Fig. 18, Buchner, 1981 (in Ali)
delta_phi = 4.6*pi/180; % angle between adjacent ommatidia
delta_rho = delta_phi*1.1; %
theta = -pi:pi/480:pi - pi/480;

% From Snyder (1979) as cited in Burton & Laughlin (2003)
filt =  exp( -4.*log(2).*abs(theta).^2 ./ delta_rho^2 ); % spatial blurring filter

eye_filt(:,33) = filt;
cnt = 1;
for j = 34:64 % simulate 64 ommatidia
    eye_filt(:,j) = circshift(filt, [0 cnt*15]);
    cnt = cnt + 1;
end
cnt = 1;
for j = 32:-1:1 % simulate 64 ommatidia
    eye_filt(:,j) = circshift(filt, [0 -cnt*15]);
    cnt = cnt + 1;
end

eye_filt(find(eye_filt < 0.005)) = 0;  % set vey low values to
% zero, so may be able to use eye_filt as a sparse matrix...not a time
% saver, because gets mulitples with non-sparse matrices.

[num_samp_pts, num_receptors] = size(eye_filt);

end