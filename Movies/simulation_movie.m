function [] = simulation_movie(bodyAng,headAng,refAng,tout,varargin)
%% simulation_movie:
%   INPUT:
%       bodyAng     : body angles
%       headAng     : head angles
%       refAng      : reference angles
%       tout        : time vector
%       rootdir     : directory to save movie >>> varargin{1} 
%   OUTPUT:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
% Example Input %
%---------------------------------------------------------------------------------------------------------------------------------
% clear ; clc ; close all
% tout = 0:(1/100):5;
% bodyAng = 40*sin(2*pi*1*tt + 0);
% headAng = 40*sin(2*pi*2*tt + 0);
% rootdir = 'Q:\Box Sync\Research';
%---------------------------------------------------------------------------------------------------------------------------------
vidFs = 100; % video frame rate [Hz]
if nargin==5
    export = true;
    root.mov = [varargin{1} '\mov'];
    [status,msg,~] = mkdir(root.mov); % create directory for .mat files
    if status
        warning(msg)
        disp(['Folder located: ' root.mov])
    else
        error('Directory not created')
    end
    VID = VideoWriter([root.mov '\TEST' '.avi'],'Uncompressed AVI'); % create video object
    VID.FrameRate = vidFs;
    open(VID)
elseif nargin<5
    export = false;
else
    error('Too may inputs')
end

% Body geometry
body.center = [0,0];
body.L = 40;
body.ratio = 0.7;
body.ecc = 0.9;
body.C = [0 0 1];

% Head geometry
head.L = 15;
head.ratio = 0.5;
head.ecc = 0.2;
head.C = [1 0 0];

% Time for video
Ts = 1/vidFs;
tt = tout(1):(1/vidFs):tout(end);

% Transform kinematics into video time
body.ang = interp1(tout, bodyAng, tt, 'nearest'); % interpolate body for video
head_ang = interp1(tout, headAng, tt, 'nearest'); % interpolate head for video
head.ang = body.ang + head_ang; % shift head angle to global coordinate frame
ref.ang  = interp1(tout, refAng, tt, 'nearest'); % interpolate pattern for video
ref.pos  = round((96/360)*ref.ang);
ref.pos(ref.pos<=0) = ref.pos(ref.pos<=0) + 96;

FIG = figure (1);
FIG.Color = 'k';
FIG.Position = [400 50 900 900];
subplot(12,1,1:8) ; cla ; hold on ; axis square ; axis equal ; axis off ; axis(60*[-1 1 -1 1])
subplot(12,1,9:12) ; cla ; hold on
plot([0 tt(end)],[0 0],'Color',[0.5 0.5 0.5],'LineStyle','--')
h1 = animatedline('Color','r','LineWidth',2); % for head angle
h2 = animatedline('Color','b','LineWidth',2); % for body angle
h3 = animatedline('Color',[0.5 0.2 0.9],'LineWidth',2,'LineStyle','-'); % for gaze angle
h4 = animatedline('Color','g','LineWidth',2,'LineStyle','--'); % for pattern angle

tic
for kk = 1:length(body.ang)
    subplot(12,1,1:8) ; hold on
%     [~] = draw_pattern(ref.pos(kk),7,pattern,body.center,50,5);    
	[h.body,body.top,~] = draw_ellipse(body.center,body.L,body.ratio,body.ecc,body.ang(kk),body.C);
    [h.head,~,~] = draw_semi_ellipse(body.top,head.L,head.ratio,head.ecc,head.ang(kk),head.C);
    [h.ref,~,~] = draw_ellipse_pattern(body.center,100,0.5,0,15,ref.ang(kk));
    
    subplot(12,1,9:12) ; hold on ; set(gca, 'color', 'w')
 	ylabel('deg','Color','w','FontSize',15)
	xlabel('time (s)','Color','w','FontSize',15);
    xlim([0 round(tt(end))])
    maxY = 10*round(max(abs([ref.ang,head.ang]))/10);
    intY = maxY/5;
    ylim(1.1*maxY*[-1 1])
 	set(gca,'ycolor','w');
    set(gca,'xcolor','w');
    set(gca,'YTick',(-maxY:intY:maxY))
    set(gca,'XTick',0:1:round(tt(end)))
 	addpoints(h1,tt(kk),head_ang(kk))
 	addpoints(h2,tt(kk),body.ang(kk))
 	addpoints(h3,tt(kk),head.ang(kk))
   	addpoints(h4,tt(kk),ref.ang(kk))
    
    pause(Ts/2)
    
    if export
        writeVideo(VID,getframe(FIG)); % write frame to .avi
    end
    
    if kk~=length(body.ang)
        for jj = 1:length(h.body)
            delete(h.body{jj})
        end
        for jj = 1:length(h.head)
            delete(h.head{jj})
        end
        for jj = 1:length(h.ref)
            delete(h.ref{jj})
        end
    end
end
toc

disp('DONE')
end