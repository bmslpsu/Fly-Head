function [EMD_ALL, Pattern_ALL, Eye_ALL] = EMDsim_v5_7_5deg_Grnd_new(showplot, export)
% This simulation reconstructs the fly's visual environment and calculates
% optic flow based on Reichart detector

%close all;

% % Folder structure
% dirXY = [root 'x_y\'];
% dirCent = [root 'centroid\'];
% dirThresh = [root 'threshold\'];
% dirVid = [root 'vid\'];
% dirAng = [root 'angles\'];
% dirPos = [root 'pos\'];

dirPat = 'D:\JMM\MatlabRoot_Magno\Panels\Patterns\';
% Open pattern file
[files, dirPat] = uigetfile({'*.mat', 'MAT-files'}, ...
    'Select corresponding pattern file', dirPat, 'MultiSelect','off');
%
% open file
patfile = [dirPat files];
load(patfile);

% make directory for video
mkdir('C:/JMM/Magno_data/7_5deg_Grnd/', 'sim5');

% make eye filters
[eye_filt, num_samp_pts, num_receptors] = make_eye_filters();

% TEMPORALLY FILTER SPACE-TIME EYE SIGNAL
t  = 0:0.001:0.060;
tp = 0.020; % 20 ms time to peak; from Dickson et al. 2008
sigma = 0.355; % width
% linear impulse response light energy -> neural signal
Vt = exp(-((log(t./tp).^2) ./ (2*sigma^2)));

% initializations for HR model
samp_rate = 100;
tau = 0.035;
h = 1/samp_rate;
% discrete low-pass filter parameters
a = h / (tau+h);

InMat = 5*(rand(1,num_receptors) - 0.5);
FiltMat = zeros(size(InMat)); FiltMat_1 = zeros(size(InMat));

% specify position
step = 10; % interval for 37deg angular vel for sampling of 100
Pos = ones(step*96,1); % initialize
tt = 0:1/samp_rate:length(Pos)/samp_rate-1/samp_rate; % time index
Pos(1:100) = 1; % init start at 1
for jj = 1:96-1
    
    Pos(jj*step:(jj+1)*step) = jj;
    
end

% preallocate array
EMD_ALL = nan(length(Pos), num_receptors);
Pattern_ALL = nan(length(Pos), 96); % 96-pixel column
Eye_ALL = nan(length(Pos), num_receptors);

mm = 1;
for jj = 1:length(Pos)
    
     if showplot == 1
        
        H = figure(50); cla; clf;
        set(H, 'Renderer','OpenGL');
        set(H, 'Position',[100, 100, 16*40, 16*50]);
        %set(H, 'Position',[100, 100, 16*20, 16*30]);
        
     end
    
     pat = pattern.Pats(:, :, Pos(jj));
     
     if showplot == 1
         x = [-180:3.75:-3.75 3.75:3.75:180];
         y = [-30:3.75:-3.75 3.75:3.75:30];
         hold on;
         mymap = [zeros(64,1) linspace(0, 1, 64)' zeros(64,1)]; % custom color map
         
         subplot(4,1,1:2)
         imagesc(x, y, pat) % show display
         %subimage(x, y, pat, mymap) % show display
         colormap(mymap);
         xlim([-180 180]); ylim([-30 30]);
         set(gca,'XTick',[-180 0 180]); set(gca,'YTick',[-30 0 30]);
         set(gca, 'YTickLabel', {'30','0', '-30'})
         set(gcf, 'color', 'k');
         set(gca, 'color', 'k')
         set(gca,'xcolor','w','ycolor','w');
         set(gca,'box','off');
         freezeColors;
         
     end   
     
     
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATE OUTPUT FROM ELEMENTARY MOTION DETECTOR
    % (HASSENTSEIN-REICHARDT)
    
    current_frame = squeeze(pattern.Pats(:, :, Pos(jj)));
    % upsample by factor of 10
    for k = 1:10
        Up_frame(k:10:num_samp_pts) = current_frame(1,1:96);
    end
    
    % get filtered eye projection
    eye_sample = Up_frame*eye_filt;    
    
    % compute HR motion
    InMat = eye_sample;
    FiltMat = a*(InMat) + (1-a)*FiltMat_1; % discrete low-pass filter
    FiltMat_1 = FiltMat;
    % HR_Motion = Va(t-tau) * Vb - Vb(t-tau) * Va
    HR_Motion = (FiltMat(1:end-1).*InMat(2:end) - FiltMat(2:end).*InMat(1:end-1)); % delay and correlate
    
    if showplot == 1
        
        subplot(4,1,3:4)
        plot(HR_Motion,'w');
        xlim([1 64]); ylim([-300 300]);
        set(gca,'XTick',[1 64]); set(gca,'YTick',[-300 0 300]);
        set(gca, 'YTickLabel', {'','', ''})
        set(gca, 'XTickLabel', {'1','64'})
        set(gcf, 'color', 'k');
        set(gca, 'color', 'k')
        set(gca,'xcolor','w','ycolor','w');
        xlabel('ommatidia')
        ylabel('EMD output')
        set(gca,'box','off');
        % display time stamp
        text(30, -430, sprintf('%0.3f',tt(jj)),'fontsize',10,'color',[1 1 1])
        
    end
    
    % prepare data for space-time plot
    EMD_ALL(jj,1:end-1) = HR_Motion;
    Pattern_ALL(jj,:) = current_frame(1,1:96);
    Eye_ALL(jj, :) = eye_sample;
    
    % get frame for movie
    %mov(mm) = getframe;
    
    % export image
    if export == 1
        imgdirUnix = ['C:/JMM/Magno_data/7_5deg_Grnd/sim5/'];
        filename = sprintf([imgdirUnix 'image%04d.jpg'], mm);
        export_fig(gcf, filename, '-q95','-nocrop');
        
    end
    
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