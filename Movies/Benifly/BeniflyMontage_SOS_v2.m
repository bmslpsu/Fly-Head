function [MOV] = BeniflyMontage_SOS_v2(rootdir,rootpat,vidFs,export)
%% BeniflyMontage_SOS_v2:  makes movie for fly in rigid tether
%                   Includes fly video, head tracking, wing tracking,pattern position and plots of data
%   INPUT:
%       rootdir     : directory containing BENIFLY file
%       rootpat     : directory containing PATTERN file
%       vidFs       : video display FPS
%       export      : boolean (1=export video to images)
%   OUTPUT:
%       MOV         : structure containing movie 
%
% Example Input %
% clear ; clc ; close all
% export = true;
% vidFs = 50;
% rootdir = 'H:\EXPERIMENTS\Experiment_SOS_v2\Benifly\new';
% rootdir = 'H:\EXPERIMENTS\Experiment_Sinusoid\15\Vid\Benifly\new';
% rootpat = 'C:\Users\boc5244\Documents\GitHub\Arena\Patterns';

% Select angle file
[FILE.benifly, PATH.benifly] = uigetfile({'*.csv', 'DAQ-files'}, ...
    'Select ANGLE file', rootdir, 'MultiSelect','off');

% Select pattern file
[FILE.pat, ~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select PATTERN file', rootpat, 'MultiSelect','off');

% Set file names
FILE.basename = FILE.benifly(1:end-4);
FILE.daq = [FILE.basename '.mat'];
FILE.vid = [FILE.basename '.avi'];
FILE.montage = [FILE.basename '_Montage_v2.avi'];

% Make path for pattern positions from daq and video time (assume one folder back)
pathparts = strsplit(PATH.benifly,filesep);
PATH.daq = fullfile(pathparts{1:end-4});
PATH.vid = fullfile(pathparts{1:end-3});

% Load data
disp('Loading Data...')
pattern_data = load(fullfile(rootpat,FILE.pat),'pattern'); % load pattern
benifly_data = ImportBenifly(fullfile(rootdir,FILE.benifly)); % load Benifly tracked kinematics
daq_data = load(fullfile(PATH.daq,FILE.daq),'data','t_p'); % load daq oattern positions
vid_data = load(fullfile(PATH.vid,FILE.daq),'t_v'); % load raw video
vidRead = VideoReader(fullfile(rootdir,FILE.vid));
benifly_vid = read(vidRead,[1 Inf]);
disp('DONE')

% Create directories
root.mov = [PATH.benifly '\Movie']; % movie directory
mkdir(root.mov) % create directory for export images
% root.image = [root.daq 'Movie\' dirName]; % image directory

% Get video, pattern, position, & angles data
Fly.vid = benifly_vid; % video data
Fly.time = vid_data.t_v; % video time
Fly.Fs = round(1/mean(diff(Fly.time))); % video sampling rate
Fly.Fc = 15; % cut off frequency for lpf
[b,a] = butter(2,Fly.Fc/(Fly.Fs/2),'low'); % make lpf
Fly.head = filtfilt(b,a,rad2deg(benifly_data.Head));
Fly.wba = rad2deg(benifly_data.LWing - benifly_data.RWing);
Fly.wba = hampel(Fly.time,Fly.wba);
Fly.wba = filtfilt(b,a,Fly.wba - mean(Fly.wba));

[Fly.xP,Fly.yP,Fly.dp,n_frame] = size(Fly.vid ); % get size of video
center = [round(Fly.yP/2) , round(Fly.xP/2)]; % center point for pattern & fly
radius.center = floor(max([Fly.yP Fly.xP])/1.45); % radius of pattern
radius.width = 12; % radius display width
rin  = radius.center - radius.width;
rout = radius.center + radius.width;
x1 = center(1);
y1 = center(2);
sA = 3.75 * pi/180; % angle pixel subtends
Pat.pos = round((96/10)*(daq_data.data(:,2)-mean(0))); % pattern position
Pat.time = daq_data.t_p; % pattern time
Pat.int = interp1(Pat.time, Pat.pos, Fly.time, 'nearest'); % interpolate pattern to match fly video
pat_lim = 40;
Pat.wrap = wrapdata(Pat.int,pat_lim,false);
Stim = 3.75*(Pat.wrap - mean(Pat.wrap));

% Create structure to store frames
MOV(1:n_frame) = struct('cdata', [], 'colormap',[]);

% Create video object
if export
    VID = VideoWriter(fullfile(root.mov,FILE.montage),'Uncompressed AVI');
    VID.FrameRate = vidFs;
    open(VID)
end

FIG = figure (1); clf % main figure window for display & export
set(FIG, 'color', 'k');
set(FIG, 'Renderer','OpenGL');
set(FIG, 'Position',0.8*[100, 100, 16*40, 16*50]);
clear ax h

subplot(14,1,1:8) ; cla ; hold on; axis square % for fly & pattern vid
linewidth = 0.75;
subplot(14,1,9:10)  ; cla ; hold on ; h.pat = animatedline('Color','g','LineWidth',linewidth); % for pattern angle
subplot(14,1,11:12) ; cla ; hold on ; h.head = animatedline('Color','b','LineWidth',linewidth); % for head angle
subplot(14,1,13:14) ; cla ; hold on ; h.wba = animatedline('Color','r','LineWidth',linewidth); % for wba angle

pp = 1;
iter = round(Fly.Fs/vidFs);
disp('Exporting Video...')
for jj = 1:n_frame % for each frame
    expframe = mod(jj,iter)==0;
    if expframe % if we want to display this frame
        pat = pattern_data.pattern.Pats(1,:,round(Pat.int(jj)),5); % top row of pattern
        patS = circshift(pat,[0 0]); % shift pattern to fly reference frame

        I = find(patS~=0);
        theta = (I.*3.75) .* (2*pi/360); % lit pixels
        theta_ALL = deg2rad(3.75*(1:96));

        Frame = Fly.vid(:,:,:,jj); % current raw frame
        DISP = Frame(40:end,:,:); % video frame to display

        % Display fly video
        subplot(14,1,1:8) ; cla ; hold on; axis square
        imshow(DISP); hold on

        % Make pattern ring
        for kk = 1:length(theta_ALL)
            xin = x1 + rin*cos(theta_ALL(kk));
            xout = x1 + rout*cos(theta_ALL(kk));
            xinN = x1 + rin*cos(theta_ALL(kk) + sA);
            xoutN = x1 + rout*cos(theta_ALL(kk) + sA);
            yin = y1 + rin*sin(theta_ALL(kk));
            yout = y1 + rout*sin(theta_ALL(kk));
            yinN = y1 + rin*sin(theta_ALL(kk) + sA);
            youtN = y1 + rout*sin(theta_ALL(kk) + sA);

            if sum(ismember(theta, theta_ALL(kk))) == 1 % if lit
                patch([xout, xoutN, xinN, xin], [yout, youtN,yinN, yin],'g','linestyle','none',...
                    'FaceAlpha',pat(kk)*(1/(2^(pattern_data.pattern.gs_val)-1)));
            else % if dark
                patch([xout, xoutN, xinN, xin], [yout, youtN,yinN, yin],'k','linestyle','none');
            end
        end
    end
    
    % Pattern plot
 	ax(1) = subplot(14,1,9:10) ; hold on
    fsize = 10;
 	ylabel('Stimulus (°)','Color','w','FontSize',fsize);
    addpoints(h.pat, Fly.time(jj), Stim(jj))
    
    % Head plot
    ax(2) = subplot(14,1,11:12) ; hold on
	ylabel('Head (°)','Color','w','FontSize',fsize)
    addpoints(h.head, Fly.time(jj), Fly.head(jj))
    
   	% Wing plot
    ax(3) = subplot(14,1,13:14) ; hold on
	ylabel('\Delta WBA (°)','Color','w','FontSize',fsize)
  	xlabel('Time (s)','Color','w','FontSize',fsize)    
    addpoints(h.wba, Fly.time(jj), Fly.wba(jj))
    
    drawnow
    
    set(ax, 'Color', 'k', 'YColor', 'w', 'XColor', 'w', 'XLim', [0 round(Fly.time(end))]) 
    set(ax(1:2), 'YLim', 20*[-1 1], 'YTick', [-15 0 15])
    set(ax(3), 'YLim', 30*[-1 1], 'YTick', [-25 0 25])
    set(ax(1:2), 'XTick', [])
    set(ax(3),  'XTick', 0:2:round(Fly.time(end)))
    set(ax, 'FontSize', 9)
    
    if expframe
        % Store frame
        MOV(pp) = getframe(FIG);

        if export
            writeVideo(VID,getframe(FIG));
        end
        pp = pp + 1;
    end
end
disp('DONE')
disp('Saving...')
if export
    close(VID) % close .avi
    Fs = Fly.Fs;
% 	save([root.mov dirName '.mat'],'MOV','Fs','-v7.3','-nocompression') % save movie as .mat file
end
disp('DONE')
end