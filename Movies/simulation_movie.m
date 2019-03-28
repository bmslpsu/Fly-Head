function [MOV] = simulation_movie(rootdir,rootpat,export,vidFs)
%% simulation_movie: on
%   INPUT:
%       rootdir     : export directory
%       rootpat     : directory containing PATTERN files
%       vidFs       : video display FPS
%       export      : boolean (1=export video to images)
%   OUTPUT:
%       MOV         : structure containing movie 
%---------------------------------------------------------------------------------------------------------------------------------
% Example Input %
% clear ; clc ; close all
% export = true;
% vidFs = 50;
% rootpat = 'Q:\Box Sync\Git\Arena\Patterns\';
%---------------------------------------------------------------------------------------------------------------------------------
% Set directories
root.pat    =  rootpat; % pattern location

% Select pattern file
[FILE.pat, ~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select PATTERN file', root.pat, 'MultiSelect','off');

% Load data
disp('Loading Data...')
data = [];
load([root.pat  FILE.pat],'pattern') % load pattern
disp('DONE')

% Create directories
[~,dirName,~] = fileparts([rootdir 'TEST_v0']); % get file name
root.mov = [rootdir 'Movie\']; % movie directory
root.image = [rootdir 'Movie\' dirName]; % image directory
mkdir(root.image) % create directory for export images
%%
close all;clc

body.e = 0.9;
head.e = 0.2;

body.L = 10;
head.L = 4;

r = 0.7;
body.center = [0,0];
ang.body = (0:360);
% tt = 0:0.01:10;
% ang.body = 45*sin(2*pi*1*tt);
ang.head = 0*linspace(ang.body(1),ang.body(end),length(ang.body));

figure (1) ; clf ; hold on ; axis equal
set(gcf,'Color','w')
axis off
axis(10*[-1 1 -1 1])
for kk = 1:length(ang.body)
    body.top(1) = body.center(1) + body.L*(1-r)*sind(ang.body(kk));
    body.top(2) = body.center(2) + body.L*(1-r)*cosd(ang.body(kk));
    
 	body.bot(1) = body.center(1) - body.L*r*sind(ang.body(kk));
	body.bot(2) = body.center(2) - body.L*r*cosd(ang.body(kk));
    
    xx = [body.center(1) body.top(1) body.bot(1)];
    yy = [body.center(2) body.top(2) body.bot(2)];
    
    a = 1/2*sqrt((body.top(1)-body.bot(1))^2+(body.top(2)-body.bot(2))^2);
    b = a*sqrt(1-body.e^2);
    t = linspace(0,2*pi);
    X = a*cos(t);
    Y = b*sin(t);
    w = atan2(body.bot(2)-body.top(2),body.bot(1)-body.top(1));
    x = (body.top(1)+body.bot(1))/2 + X*cos(w) - Y*sin(w);
    y = (body.top(2)+body.bot(2))/2 + X*sin(w) + Y*cos(w);
    h3 = plot(x,y,'k-','LineWidth',2);
 	h7 = patch(x,y,'b');
    
    h2 = plot(xx,yy,'-k','LineWidth',2);  
    
	head.top(1) = body.top(1) + head.L*sind(ang.head(kk));
    head.top(2) = body.top(2) + head.L*cosd(ang.head(kk));
    
	xx = [body.top(1) , body.top(1) + head.L*sind(ang.head(kk))/2];
    yy = [body.top(2) , body.top(2) + head.L*cosd(ang.head(kk))/2];
    
	a = 1/2*sqrt((head.top(1)-body.top(1))^2+(head.top(2)-body.top(2))^2);
    b = a*sqrt(1-head.e^2);
    t = linspace(pi/2,(3/2)*pi);
    X = a*cos(t);
    Y = b*sin(t);
    w = atan2(body.top(2)-head.top(2),body.top(1)-head.top(1));
    x = (body.top(1)) + X*cos(w) - Y*sin(w);
    y = (body.top(2)) + X*sin(w) + Y*cos(w);
    
    h6 = patch([x x(1)],[y y(1)],'r');
	h5 = plot([x x(1)],[y y(1)],'-k','LineWidth',2);
 	h4 = plot(xx,yy,'-k','LineWidth',2);

	h1 = plot(body.top(1),body.top(2),'-ko','MarkerSize',5,'MarkerFaceColor','k');

    h9 = plot(head.top(1),head.top(2),'-ko','MarkerSize',5,'MarkerFaceColor','g');
    
    h8 = plot(body.center(1),body.center(2),'-ko','MarkerSize',5,'MarkerFaceColor','k');

    pause(0.01)
    
	delete(h1)
    delete(h2)
    delete(h3)
  	delete(h4)
    delete(h5)
    delete(h6)
	delete(h7)
    delete(h8)
    delete(h9)
end
close














%%
% Get video, pattern, position, & angles data 
[Fly.xP,Fly.yP,nFrame] = size(Fly.vid ); % get size of video
center = [round(Fly.yP/2) , round(Fly.xP/2)+45]; % center point for pattern & fly
radius.center = floor(max([Fly.yP Fly.xP])/1.45); % radius of pattern
radius.width = 10; % radius display width
rin  = radius.center - radius.width;
rout = radius.center + radius.width;
x1 = center(1);
y1 = center(2);
sA = 3.75 * pi/180; % angle pixel subtends
Pat.pos = round((96/10)*(data(:,2)-mean(0))); % pattern position
Pat.time = t_p; % pattern time
Pat.int = interp1(Pat.time, Pat.pos, Fly.time, 'nearest'); % interpolate pattern to match fly video

% Get head rotation point if not found
if exist('hCenter','var')~=1 
    figure ; clf ; hold ; title('Select head rotation point, press space when done')
    imshow(Fly.vid(:,:,1))
    hp = impoint();
    pause
    hCenter = getPosition(hp);
end

% Create structure to store frames
MOV(1:nFrame) = struct('cdata', [], 'colormap',[]);

% Create video object
VID = VideoWriter([root.mov dirName '.avi'],'Uncompressed AVI');
VID.FrameRate = vidFs;
open(VID)

FIG = figure ; clf % main figure window for display & export
set(gcf, 'color', 'k');
set(FIG, 'Renderer','OpenGL');
set(FIG, 'Position',[100, 100, 16*40, 16*50]);
subplot(12,1,1:8) ; cla ; hold on; axis square % for fly & pattern vid
subplot(12,1,9:10)  ; cla ; hold on ; h1 = animatedline('Color','g','LineWidth',2); % for pattern angle
subplot(12,1,11:12) ; cla ; hold on ; h2 = animatedline('Color','b','LineWidth',2); % for head angle
pp = 1;
iter = round(Fly.Fs/vidFs);
disp('Exporting Video...')
for jj = 1:iter:nFrame % for each frame    
	pat = pattern.Pats(1,:,round(Pat.int(jj)),4); % body.top row of pattern
	patS = circshift(pat,[0 0]); % shift pattern to fly reference frame
    
    I = find(patS~=0);
    theta = (I.*3.75) .* (2*pi/360); % lit pixels
    theta_ALL = deg2rad(3.75*(1:96));

	Frame = Fly.vid(:,:,jj); % current raw frame
    DISP = Frame; % video frame to display
    
    % Display fly video
    subplot(12,1,1:8) ; cla ; hold on; axis square
    imshow(DISP); hold on
    hTipX = hCenter(1) + rout*sind(hAngles(jj));
    hTipY = hCenter(2) - rout*cosd(hAngles(jj));
    plot([hCenter(1),hTipX],[hCenter(2),hTipY],'-b','LineWidth',2)
    plot(hCenter(1),hCenter(2),'oc','MarkerSize',2)
    plot(x1,y1,'r.','MarkerSize',1) % display center point

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
                'FaceAlpha',pat(kk)*(1/(2^(pattern.gs_val)-1)));
        else % if dark
            patch([xout, xoutN, xinN, xin], [yout, youtN,yinN, yin],'k','linestyle','none');
        end
    end
    
    % Pattern plot
 	subplot(12,1,9:10) ; hold on ; set(gca, 'color', 'w')
 	ylabel('Display ($^{\circ}$)','Interpreter','latex','Color','w','FontSize',12);
    xlim([0 round(Fly.time(end))])
    ylim([-20 20])
 	set(gca,'ycolor','w');
    set(gca,'xcolor','k');
    set(gca,'YTick',[-20 0 20])
    set(gca,'XTick',0:1:round(Fly.time(end)))
    addpoints(h1,t_v(jj),3.75*Pat.int(jj) - 3.75*mean(Pat.int))
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
    addpoints(h2,t_v(jj),hAngles(jj) - mean(hAngles))
    drawnow
    
    % Store frame
    MOV(pp) = getframe(FIG);
    
    if export
        % Export frame to image
        filename = sprintf('image%04d.jpg', pp);
        export_fig(gcf, [root.image '\' filename], '-q95','-nocrop');
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
	save([root.mov dirName '.mat'],'MOV','Fs','-v7.3','-nocompression') % save movie as .mat file
end
disp('DONE')
end