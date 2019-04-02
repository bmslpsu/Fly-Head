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
% bodyAng = 40*sin(2*pi*1*tout + 97);
% headAng = 180*tout;
% refAng = 0*bodyAng;
% rootdir = 'C:\Users\boc5244\Box Sync\Research';
%---------------------------------------------------------------------------------------------------------------------------------
vidFs = 50; % video frame rate [Hz]
if nargin==6
    export = true;
    root.mov = varargin{1};
    filename = varargin{2};
    VID = VideoWriter(fullfile(root.mov,[filename '.avi']),'Uncompressed AVI'); % create video object
    VID.FrameRate = vidFs;
    open(VID)
elseif nargin==5
    error('Must input filename')
elseif nargin<5
    export = false;
else
    error('Too may inputs')
end

% Body geometry
body.center = [0,0];
body.L = 40;
body.ratio = 1;
body.ecc = 0.9;
body.C = [0 0 1];

% Head geometry
head.L = 0.1875*body.L;
head.ratio = 0.5;
head.ecc = 0.2;
head.C = [1 0 0];

% Stimulus geometry
ref.radius = 65;
fig.radius = 10;

% Time for video
Ts = 1/vidFs;
tt = tout(1):(1/vidFs):tout(end);

% Transform kinematics into video time
body.ang = interp1(tout, bodyAng, tt, 'nearest')'; % interpolate body for video
head.ang = interp1(tout, headAng, tt, 'nearest')'; % interpolate head for video
gaze.ang = body.ang + head.ang; % shift head angle to global coordinate frame
ref.ang  = interp1(tout, refAng, tt, 'nearest')'; % interpolate pattern for video
head.top = body.L*body.ratio*[sind(body.ang) , cosd(body.ang)];
fig.pos  = (ref.radius)*[sind(ref.ang)  , cosd(ref.ang) ];
fig.ang  = atan2d((fig.pos(:,2)-head.top(:,2)),(fig.pos(:,1)-head.top(:,1)));
fig.L    = sqrt((fig.pos(:,2)-head.top(:,2)).^(2) + (fig.pos(:,1)-head.top(:,1)).^(2));
fig.pos  = head.top + fig.L.*[cosd(fig.ang)  , sind(fig.ang)];


FIG = figure (1);
FIG.Color = 'k';
FIG.Position = [400 50 900 900];
movegui(FIG,'center')
subplot(12,1,1:8) ; cla ; hold on ; axis square ; axis equal ; axis off ; set(gca, 'color', 'w')
axis(1.4*ref.radius*[-1 1 -1 1])
subplot(12,1,9:12) ; cla ; hold on
plot([0 tt(end)],[0 0],'Color',[0.5 0.5 0.5],'LineStyle','--')
h1 = animatedline('Color','r','LineWidth',2); % for head angle
h2 = animatedline('Color','b','LineWidth',2); % for body angle
h3 = animatedline('Color',[0.5 0.2 0.9],'LineWidth',2,'LineStyle','-'); % for gaze angle
h4 = animatedline('Color','g','LineWidth',2,'LineStyle','--'); % for pattern angle

tic
for kk = 1:length(body.ang)
    subplot(12,1,1:8) ; hold on
        [h.body,body.top,~] = draw_ellipse(body.center,body.L,body.ratio,body.ecc,body.ang(kk),body.C);
        [h.head,head.top,~] = draw_semi_ellipse(body.top,head.L,head.ratio,head.ecc,gaze.ang(kk),head.C);
        [h.ref,~,~] = draw_ellipse_pattern(body.center,2*ref.radius,0.5,0,15,ref.ang(kk));
        [h.fig,~,~] = draw_ellipse(fig.pos(kk,:),fig.radius,0.5,0,ref.ang(kk),[0 0.9 0.9]);
        [h.gaze,~]  = draw_gaze(head.top,ref.radius,30,gaze.ang(kk),[0.6 0.1 0.9]);
    subplot(12,1,9:12) ; hold on
        ylabel('deg','Color','w','FontSize',15)
        xlabel('time (s)','Color','w','FontSize',15);
        xlim([0 round(tt(end))])
        maxY = 10*round(max(abs([ref.ang;gaze.ang]))/10);
        intY = maxY/5;
        ylim(1.1*maxY*[-1 1])
        set(gca,'ycolor','w');
        set(gca,'xcolor','w');
        set(gca,'YTick',(-maxY:intY:maxY))
        set(gca,'XTick',0:1:round(tt(end)))
        addpoints(h1,tt(kk),head.ang(kk))
        addpoints(h2,tt(kk),body.ang(kk))
        addpoints(h3,tt(kk),gaze.ang(kk))
        addpoints(h4,tt(kk),ref.ang(kk))
    
    pause(Ts/4)
    
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
        for jj = 1:length(h.fig)
            delete(h.fig{jj})
        end
        for jj = 1:length(h.gaze)
            delete(h.gaze{jj})
        end
    end
end
toc

disp('DONE')
end