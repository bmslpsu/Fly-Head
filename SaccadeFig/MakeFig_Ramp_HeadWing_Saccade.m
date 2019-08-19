function [FIG] = MakeFig_Ramp_HeadWing_Saccade()
%% MakeFig_Ramp_HeadWing_Saccade:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
% root = 'H:\DATA\Rigid_Data\';
root = 'Q:\Box Sync\';
[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

% load(fullfile(root,FILE{1}),'SACCADE','INTERVAL','SACD','Stim','U','N','I','TRIAL','FLY','GRAND');
load(fullfile(root,FILE{1}),'SACD','U','N','TRIAL');
%%
clc
clearvars -except SACCADE INTERVAL SACD Stim U N I TRIAL FLY GRAND

CC = repmat({[1 0 0],[0 1 0],[0 0 1]},1,2);

Vel = 3.75*U{1,3}{1};

headIdx = 2;
wingIdx = 3;

clear FIG ax
FIG = figure (1) ; clf ; hold on
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 12 4];
movegui(FIG,'center')
FIG.Name = 'Head vs Wing';
FIG.Visible = 'off';

ax.L = gca ; cla ; hold on
ax.L.FontSize = 12;
ax.L.Color = 'none';
ax.L.YColor = [0 0 0.7];
xlabel('Time (s)','FontSize',12,'FontWeight','bold')
ylabel(['Head (' char(176) '/s)'],'FontSize',12,'FontWeight','bold')
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

PP = nan(size(SACD.Head,1),3); pp = 1;
debug = false;
FIG.Visible = 'off';
for jj = 1:N.vel
    for kk = 1:N.fly
        for ii = 1:size(TRIAL{kk,jj},1)
         	HEAD = TRIAL{kk,jj}{ii,headIdx};
            WING = TRIAL{kk,jj}{ii,wingIdx};
            
          	[Head.SACD,Head.thresh,Head.count] = SacdDetect(HEAD.X(:,1),HEAD.Time,2.5,debug);
            Head.svel = HEAD.X(:,2);
            Head.svel( (Head.svel<Head.thresh)  & (Head.svel>0))  = 0;
            Head.svel( (Head.svel>-Head.thresh) & (Head.svel<0))  = 0;
            
         	[Wing.SACD,Wing.thresh,Wing.count] = SacdDetect(WING.X(:,1),WING.Time,2.5,debug);
            Wing.svel = WING.X(:,2);
            Wing.svel( (Wing.svel<Wing.thresh)  & (Wing.svel>0) )  = 0;
            Wing.svel( (Wing.svel>-Wing.thresh) & (Wing.svel<0) )  = 0;
            
            axes(ax.L) ; cla ; hold on
            plot(HEAD.Time, Head.svel,'Color', ax.L.YColor)
            plot(Head.SACD.PeakTime, Head.SACD.PeakVel, '*', 'Color', 'b', 'MarkerSize', 15)

            axes(ax.R) ; cla ; hold on
            plot(WING.Time, Wing.svel,'Color', ax.R.YColor)
            plot(Wing.SACD.PeakTime, Wing.SACD.PeakVel, '*', 'Color', 'r', 'MarkerSize', 15)
            
            window = 0.1; % search wiondow [s]
            for ww = 1:Head.count
                tpeak = Head.SACD.StartTime(ww);
                tstart = Head.SACD.StartTime(ww);
                cnd = (Wing.SACD.PeakTime>=(tpeak-window)) & (Wing.SACD.PeakTime<=(tpeak+window));
                
                if sum(cnd)==0
                    PP(pp,1) = false;
                elseif sum(cnd)==1
                    PP(pp,1) = true;
                  	plot(Wing.SACD.PeakTime(cnd), Wing.SACD.PeakVel(cnd), 'o', 'Color', 'g', 'MarkerSize', 20)
                    tdiff = Wing.SACD.StartTime(cnd) - tstart;
                    match = sign(Wing.SACD.PeakVel(cnd)) == sign(Head.SACD.PeakVel(ww));
                    PP(pp,2) = tdiff;
                    PP(pp,3) = match;
                    if abs(tdiff)>0.1
                        %pause()
                    end
                elseif sum(cnd)>1
                    PP(pp,1) = false;
                    disp('More than one saccade in window')
                    pause()
                else
                    error('WUTT')
                end
                    
            pp = pp + 1;
            end
            %pause()
        end
    end
end

disp('Linked saccade ratio:')
disp(sum(PP(:,1))./size(PP,1))
disp('Median time difference (s):')
disp(nanmedian(PP(:,2)))
disp('Mathced direction saccade ratio:')
disp(nansum(PP(:,3))./size(PP,1))

clear FIG ax
FIG = figure (2) ; clf ; hold on
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 8 4];
movegui(FIG,'center')
FIG.Name = 'Head vs Wing Time Difference';
FIG.Visible = 'on';

histogram(1000*PP(:,2), 1000*(-0.15:0.005:0.15))
ax = gca;
xlabel('Time Difference (ms)','FontSize',12)




end