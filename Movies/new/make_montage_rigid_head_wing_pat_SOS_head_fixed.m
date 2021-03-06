function [MOV] = make_montage_rigid_head_wing_pat_SOS_head_fixed(rootdir,rootpat,vidFs,export)
%% make_montage_rigid_head_wing_pat_SOS_head_fixed: makes movie for fly in rigid tether
%
% 	Includes fly video, wing tracking, & pattern position
%
%   INPUT:
%       rootdir     : directory containing .mat DAQ & VIDEO files
%       rootpat     : directory containing PATTERN file
%       vidFs       : video display FPS
%       export      : boolean (1=export video to images)
%
%   OUTPUT:
%

% Example Input
clear ; clc ; close all 
export = false;
vidFs = 50;
rootdir = 'H:\EXPERIMENTS\RIGID\Experiment_SOS_v2_HeadFixed';
rootpat = 'C:\Users\BC\Box\Git\Arena\Patterns';
pat_ypos = 5;

if ~isfolder(rootdir)
    dirflag = false;
    [rootdir,mainfile,mainext] = fileparts(rootdir);
else
    dirflag = true;
end

if ~isfolder(rootpat)
	[PATH.pat,FILE.pat,patext] = fileparts(rootpat);
    FILE.pat = [FILE.pat , patext];
else
	% Select pattern file
    [FILE.pat, PATH.pat] = uigetfile({'*.mat'}, ...
        'Select pattern file', rootpat, 'MultiSelect','off');
end

% Create data paths
PATH.raw            = rootdir;
PATH.vid            = fullfile(PATH.raw,'');
PATH.beninfly_track	= fullfile(PATH.vid,'filtvid','tracked_head_wing');
PATH.mask           = fullfile(PATH.beninfly_track,'');

if dirflag
    % Select tracked angle file (use head tracked files to select)
    [FILE.raw, PATH.head_track] = uigetfile({'*.csv'}, ...
        'Select fly file', PATH.beninfly_track, 'MultiSelect','off');
else
    FILE.raw = [mainfile , mainext];
end

% Create movie output directory
PATH.mov = fullfile(PATH.raw,'movie'); % movie directory
mkdir(PATH.mov) % create directory for export images

% Set file names
[~,FILE.basename,~] = fileparts(FILE.raw);
FILE.raw        = [FILE.basename '.mat'];
FILE.benifly   	= [FILE.basename '.csv'];
FILE.montage    = [FILE.basename '_montage_active.mp4'];
FILE.mask       = [FILE.basename '.json'];

% Load data
disp('Loading Data ...')
pattern_data    = load(fullfile(PATH.pat,FILE.pat),'pattern'); % load pattern
benifly_data    = ImportBenifly(fullfile(PATH.beninfly_track,FILE.benifly)); % load Benifly tracked kinematics
raw_data        = load(fullfile(PATH.raw,FILE.raw),'data','t_p'); % load DAQ pattern positions
vid_data        = load(fullfile(PATH.vid,FILE.raw),'vidData'); % load raw video
disp('DONE')

%% Get pattern data & sync with start of visual stimulus
pattern_total_time = 20; % [s]
reg = true; % use interp times
start_idx = nan; % use first frame
add1 = true; % add 1st frame
debug = false;

[TRIG,PAT] = sync_pattern_trigger(raw_data.t_p, raw_data.data(:,2), pattern_total_time, ... 
                                raw_data.data(:,1), reg, start_idx, add1, debug);

%% Get kinematics data
FLY.time    = TRIG.time_sync; % video time
FLY.Fs      = round(1/mean(diff(FLY.time))); % video sampling rate
FLY.Fc      = 20; % cut off frequency for lpf
[b,a]       = butter(3,FLY.Fc/(FLY.Fs/2),'low'); % make lpf
FLY.lwing   = rad2deg(hampel(FLY.time,benifly_data.LWing)); % left wing angles [deg]
FLY.rwing   = rad2deg(hampel(FLY.time,benifly_data.RWing)); % right wing angles [deg]
FLY.lwing   = filtfilt(b,a,FLY.lwing); % left wing angles [deg]
FLY.rwing   = 0 + filtfilt(b,a,FLY.rwing); % right wing angles [deg]
FLY.wba     = FLY.lwing - FLY.rwing; % delta wing-beat-amplitude [deg]
[b,a]       = butter(3,15/(FLY.Fs/2),'low'); % make lpf
FLY.wba     = filtfilt(b,a,FLY.wba); % delta wing-beat-amplitude [deg]

% Normalize fly kinematics for experimental window
FLY.int_time    = TRIG.time_sync_exp; % video time
FLY.int_lwing 	= FLY.lwing(TRIG.range);
FLY.int_rwing  	= FLY.rwing(TRIG.range);
FLY.int_wba     = FLY.wba(TRIG.range);
FLY.int_rwing  	= FLY.rwing(TRIG.range);
PAT.norm        = 3.75*round(PAT.pos_exp - mean(PAT.pos_exp));

[b,a]       = butter(3, 12/(FLY.Fs/2), 'low');
PAT.norm  	= filtfilt(b, a, PAT.norm);

%% Get video data
FLY.raw = squeeze(vid_data.vidData(:,:,TRIG.range)); % raw video data
FLY.raw_crop = FLY.raw;
FLY.nframe = size(FLY.raw,3);

[FLY.raw_yP,FLY.raw_xP,~] = size(FLY.raw_crop); % get size of raw video
FLY.raw_center = [round(FLY.raw_xP/2) , 1.25*round(FLY.raw_yP/2)]; % center point for pattern & fly
radius = floor(max([FLY.raw_yP FLY.raw_xP])/1.0); % radius of pattern
thickness = 8; % radius display width

%% Get benifly parameters/mask
% Heading
fid = fopen(fullfile(PATH.mask, FILE.mask));
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
params = jsondecode(str);

FLY.body_reg = -90 - rad2deg(atan2(params.gui.head.hinge.y - params.gui.abdomen.hinge.y, ...
                            params.gui.head.hinge.x - params.gui.abdomen.hinge.x));
FLY.wing_length = 150;

FLY.lwing_hinge = [params.gui.left.hinge.x  , params.gui.left.hinge.y];
FLY.rwing_hinge = [params.gui.right.hinge.x , params.gui.right.hinge.y];

FLY.lwing_tip = FLY.lwing_hinge - FLY.wing_length*[cosd(FLY.int_lwing - FLY.body_reg),  ...
                    sind(FLY.int_lwing - FLY.body_reg)];
FLY.rwing_tip = FLY.rwing_hinge + FLY.wing_length*[cosd(FLY.int_rwing - FLY.body_reg), ...
                    -sind(FLY.int_rwing - FLY.body_reg)];

%% Make Movie
% Create structure to store frames
MOV(1:FLY.nframe) = struct('cdata', [], 'colormap',[]);

% Create video object
if export
    VID = VideoWriter(fullfile(PATH.mov,FILE.montage),'MPEG-4');
    VID.FrameRate = vidFs;
    open(VID)
end

FIG = figure (1); clf
set(FIG, 'Color', 'k', 'Renderer', 'OpenGL', 'Units', 'inches', ...
    'Position', [-16 2.5 8 10*(12/14)]);
movegui(FIG, 'center')
linewidth = 1.25;
fontsize = 12;

gs = pattern_data.pattern.gs_val + 1;
cmap = [zeros(gs,1), linspace(0,1,gs)', zeros(gs,1)];
colormap(cmap)

clear ax
ax(1) = subplot(12,1,1:7) ; cla; hold on; axis square % for raw fly & pattern vid
ax(2) = subplot(12,1,9:10) ; cla; hold on
        ylabel('Stimulus (�)','Color','w','FontSize',fontsize)
        h.pat = animatedline('Color','g','LineWidth',linewidth); % for pattern angle
ax(3) = subplot(12,1,11:12) ; cla; hold on
        ylabel('\DeltaWBA (�)','Color','w','FontSize',fontsize)
    	xlabel('Time (s)','Color','w','FontSize',fontsize)
        h.wba = animatedline('Color','r','LineWidth',linewidth); % for dWBA angle

set(ax(2:end), 'FontSize', 12, 'Color', 'k', 'YColor', 'w', 'XColor', 'w', 'FontWeight', 'bold',...
    'LineWidth', 1.5,'XLim', [-0.2 round(FLY.int_time(end))])
set(ax(end),'XTick', 0:2:round(FLY.time(end)))
set(ax(2), 'XTickLabel', [], 'XColor', 'none')
set(ax(2:end), 'YLim', 20*[-1 1], 'YTick', 15*[-1 0 1])

w_ylim = 5*ceil(max(abs(FLY.int_wba)) / 5);
set(ax(3), 'YLim', w_ylim*[-1 1], 'YTick', (w_ylim-5)*[-1 0 1])

linkaxes(ax(2:end),'x')
align_Ylabels_ax(ax(2:end)')

iter = round(FLY.Fs/vidFs); % # of frames to skip per iteration to acheive desired frame rate
expframe = circshift(mod(1:FLY.nframe,iter)==1,0); % which frames to export
disp('Exporting Video...')
tic
pat_image = 255*pattern_data.pattern.Pats(1,:,1,pat_ypos);
for jj = 1:FLY.nframe % for each frame
    if expframe(jj) % if we want to display this frame
        % Get frames
        disp(jj)
        if jj >= iter
            win = jj-(iter-1):jj;
        else
            win = jj;
        end
        Frame.raw = 1.5*median(FLY.raw(:,:,win),3); % current raw frame median across frames
        Frame.raw = medfilt2( imadjust(Frame.raw), [3 3]);
        
        % Display raw video
        subplot(12,1,1:7); cla; hold on; axis image
            imshow(Frame.raw)                         
          	plot([FLY.lwing_hinge(1) , mean(FLY.lwing_tip(win,1))], ... % update line drawn to left wing
                [FLY.lwing_hinge(2) , mean(FLY.lwing_tip(win,2))], 'Color', 'r', 'LineWidth', 2.5)
            
        	plot([FLY.rwing_hinge(1) , mean(FLY.rwing_tip(win,1))], ... % update line drawn to right wing
                [FLY.rwing_hinge(2) , mean(FLY.rwing_tip(win,2))], 'Color', 'r', 'LineWidth', 2.5)
            
            plot(FLY.lwing_hinge(1), FLY.lwing_hinge(2), 'r.', 'MarkerSize',20)
            plot(FLY.rwing_hinge(1), FLY.rwing_hinge(2), 'r.', 'MarkerSize',20)
            
            offset = [120 120];
            axis([-offset(1), FLY.raw_xP + offset(1), -offset(2), FLY.raw_yP + offset(2)])
            
            % Make pattern ring
            ax_pat = axes; axis image
            set(ax_pat, 'Color', 'none', 'XColor', 'none', 'YColor', 'none', ...
                            'Position', ax(1).Position)
           	cla

            pat_pos = 1*round(mean(PAT.norm(win)));
            theta = -deg2rad(pat_pos) + -(pi/2 + linspace(0, 2*pi, size(pat_image,2)));
            x = radius * cos(theta) + FLY.raw_center(1);
            y = radius * sin(theta) + FLY.raw_center(2);
            z = zeros(1,length(x));
            hs = surface(ax_pat,[x;x],[y;y],[z;z],[pat_image;pat_image], ...
                'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', thickness);
    end
    
   	% Stimulus plot
	subplot(12,1,9:10); hold on
        addpoints(h.pat, FLY.int_time(jj), PAT.norm(jj))
        
   	% WBA plot
	subplot(12,1,11:12); hold on
        addpoints(h.wba, FLY.int_time(jj), FLY.int_wba(jj))
        
    drawnow
    
    if export
        if expframe(jj)
            fig_frame = getframe(FIG);
         	writeVideo(VID,fig_frame);
        end
    end
    pause(0.001)
end
toc

if export
 	disp('Saving...')
    pause(1)
    close(VID) % close video
end
disp('DONE')
end