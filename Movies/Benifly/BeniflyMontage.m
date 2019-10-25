function [MOV] = BeniflyMontage(rootdir,rootpat,vidFs,export)
%% BeniflyMontage:  makes movie for fly in rigid tether
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
% rootdir = 'H:\EXPERIMENTS\Experiment_Asymmetry_Control_Verification\HighContrast\30\Vid\tracked';
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
FILE.montage = [FILE.basename '_Montage.avi'];

% Make path for pattern positions from daq and video time (assume one folder back)
pathparts = strsplit(PATH.benifly,filesep);
PATH.daq = fullfile(pathparts{1:end-3});
PATH.vid = fullfile(pathparts{1:end-2});

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
Fly.head = rad2deg(benifly_data.Head);
Fly.Fs = round(1/mean(diff(Fly.time)));
[Fly.xP,Fly.yP,Fly.dp,n_frame] = size(Fly.vid ); % get size of video
center = [round(Fly.yP/2) , round(Fly.xP/2)]; % center point for pattern & fly
radius.center = floor(max([Fly.yP Fly.xP])/1.45); % radius of pattern
radius.width = 10; % radius display width
rin  = radius.center - radius.width;
rout = radius.center + radius.width;
x1 = center(1);
y1 = center(2);
sA = 3.75 * pi/180; % angle pixel subtends
Pat.pos = round((96/10)*(daq_data.data(:,2)-mean(0))); % pattern position
Pat.time = daq_data.t_p; % pattern time
Pat.int = interp1(Pat.time, Pat.pos, Fly.time, 'nearest'); % interpolate pattern to match fly video
pat_lim = 20;
Pat.wrap = wrapdata(Pat.int,pat_lim,false);
% Pat.wrap = Pat.wrap - Pat.wrap(1);

% Create structure to store frames
MOV(1:n_frame) = struct('cdata', [], 'colormap',[]);

% Create video object
if export
    VID = VideoWriter(fullfile(root.mov,FILE.montage),'Uncompressed AVI');
    VID.FrameRate = vidFs;
    open(VID)
end

FIG = figure (1); clf % main figure window for display & export
set(gcf, 'color', 'k');
set(FIG, 'Renderer','OpenGL');
set(FIG, 'Position',0.8*[100, 100, 16*40, 16*50]);
subplot(12,1,1:8) ; cla ; hold on; axis square % for fly & pattern vid
subplot(12,1,9:10)  ; cla ; hold on ; h1 = animatedline('Color','g','LineWidth',1.5); % for pattern angle
subplot(12,1,11:12) ; cla ; hold on ; h2 = animatedline('Color','b','LineWidth',1.5); % for head angle
pp = 1;
iter = round(Fly.Fs/vidFs);
disp('Exporting Video...')
for jj = 1:iter:n_frame % for each frame 
	pat = pattern_data.pattern.Pats(1,:,round(Pat.int(jj)),5); % top row of pattern
	patS = circshift(pat,[0 0]); % shift pattern to fly reference frame
    
    I = find(patS~=0);
    theta = (I.*3.75) .* (2*pi/360); % lit pixels
    theta_ALL = deg2rad(3.75*(1:96));

	Frame = Fly.vid(:,:,:,jj); % current raw frame
    DISP = Frame; % video frame to display
    
    % Display fly video
    subplot(12,1,1:8) ; cla ; hold on; axis square
    imshow(DISP); hold on
%     hTipX = hCenter(1) + rout*sind(hAngles(jj));
%     hTipY = hCenter(2) - rout*cosd(hAngles(jj));
%     plot([hCenter(1),hTipX],[hCenter(2),hTipY],'-b','LineWidth',2)
%     plot(hCenter(1),hCenter(2),'oc','MarkerSize',2)
%     plot(x1,y1,'r.','MarkerSize',1) % display center point

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
    
    % Pattern plot
 	subplot(12,1,9:10) ; hold on ; set(gca, 'color', 'w')
 	ylabel('Display ($^{\circ}$)','Interpreter','latex','Color','w','FontSize',12);
    xlim([0 round(Fly.time(end))])
    ylim([0 3.75*pat_lim])
 	set(gca,'ycolor','w');
    set(gca,'xcolor','k');
    set(gca,'XTick',0:1:round(Fly.time(end)))
    addpoints(h1,Fly.time(jj), 3.75*Pat.wrap(jj))
    drawnow
    
    % Head plot
    subplot(12,1,11:12) ; hold on ; set(gca, 'color', 'w')
	ylabel('Head ($^{\circ}$)','Interpreter','latex','Color','w','FontSize',12)
  	xlabel('Time (s)','Interpreter','latex','Color','w','FontSize',12)    
    xlim([0 round(Fly.time(end))])
    ylim([-20 20])
 	set(gca,'ycolor','w');
 	set(gca,'xcolor','w');
    set(gca,'YTick',[-20 0 20])
    set(gca,'XTick',0:1:round(Fly.time(end)))
    addpoints(h2,Fly.time(jj),Fly.head(jj))
    drawnow
    
    % Store frame
    MOV(pp) = getframe(FIG);
    
    if export
        % Export frame to image
%         filename = sprintf('image%04d.jpg', pp);
%         export_fig(gcf, [root.image '\' filename], '-q95','-nocrop');
        % Write frame to .avi
        writeVideo(VID,getframe(FIG));
    end
    pp = pp + 1;
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