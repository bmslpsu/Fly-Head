function [FIG] = MakeFig_Ramp_Head_Saccade()
%% MakeFig_Ramp_Head_Saccade:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

load(fullfile(root,FILE{1}),'SACCADE','INTERVAL','SACD','Stim','U','N','I','TRIAL','FLY','GRAND');

CC = repmat({[1 0 0],[0 1 0],[0 0 1]},1,2);

Vel = 3.75*U{1,3}{1};

clearvars -except SACCADE INTERVAL SACD Stim U N I TRIAL FLY GRAND CC Vel
%%
clc

time = linspace(0,10,2000)';
headIdx = 2;
wingIdx = 3;

clear FIG ax
FIG = figure (1) ; clf ; hold on
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 12 4];
movegui(FIG,'center')
FIG.Name = 'Head vs Wing';

ax.L = gca ; cla ; hold on
ax.L.FontSize = 12;
ax.L.Color = 'none';
ax.L.YColor = [0 0 0.7];
xlabel('Time (s)','FontSize',12,'FontWeight','bold')
ylabel(['Head (' char(176) ')'],'FontSize',12,'FontWeight','bold')
ax.L.YLim = 1000*[-1 1];
ax.L.XLim = [0 10];

ax.R = axes;
ax.R.YAxisLocation = 'right';
ax.R.Color = 'none';
ax.R.YColor = [0.7 0 0];
ax.R.YLabel.String = 'Wing (V/s)';
ax.R.YLabel.FontSize = 12;
ax.R.YLabel.FontWeight = 'bold';
ax.R.XLim = ax.L.XLim;
ax.R.YLim = 200*[-1 1];
ax.R.XTickLabels = '';
ax.R.XColor = 'none';
ax.R.Position = ax.L.Position;
linkaxes([ax.L ax.R],'x');
% linkaxes([ax.L ax.R],'y');
debug = false;
for jj = 1:N.vel
    for kk = 1:N.fly
        for ii = 1:size(TRIAL{kk,jj},1)
         	Head.time = TRIAL{kk,jj}{ii,headIdx}.Time;
            Wing.time = TRIAL{kk,jj}{ii,wingIdx}.Time;
            Head.pos = TRIAL{kk,jj}{ii,headIdx}.X(:,1);
            Wing.pos = TRIAL{kk,jj}{ii,wingIdx}.X(:,1);
          	Head.vel = TRIAL{kk,jj}{ii,headIdx}.X(:,2);
            Wing.vel = TRIAL{kk,jj}{ii,wingIdx}.X(:,2);
            
          	[Head.SACD,Head.thresh,Head.count] = SacdDetect(Head.pos,Head.time,2.5,debug);
            Head.svel = Head.vel;
            Head.svel( (Head.vel<Head.thresh)  & (Head.vel>0))  = 0;
            Head.svel( (Head.vel>-Head.thresh) & (Head.vel<0))  = 0;
            
         	[Wing.SACD,Wing.thresh,Wing.count] = SacdDetect(Wing.pos,Wing.time,2.5,debug);
            Wing.svel = Wing.vel;
            Wing.svel( (Wing.vel<Wing.thresh)  & (Wing.vel>0))  = 0;
            Wing.svel( (Wing.vel>-Wing.thresh) & (Wing.vel<0))  = 0;
            
            axes(ax.L) ; cla ; hold on
            plot(Head.time, Head.svel,'Color', ax.L.YColor)
            plot(Head.SACD.PeakTime, Head.SACD.PeakVel, '*', 'Color', 'b', 'MarkerSize', 15)

            axes(ax.R) ; cla ; hold on
            plot(Wing.time, Wing.svel,'Color', ax.R.YColor)
            plot(Wing.SACD.PeakTime, Wing.SACD.PeakVel, '*', 'Color', 'r', 'MarkerSize', 15)

            pause()
        end
    end
end

end