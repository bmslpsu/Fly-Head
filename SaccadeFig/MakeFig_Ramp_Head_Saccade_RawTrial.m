function [FIG] = MakeFig_Ramp_Head_Saccade_RawTrial()
%% MakeFig_Ramp_Head_Saccade_RawTrial:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

load(fullfile(root,FILE{1}),'Stim','U','N','I','T','TRIAL');

Vel = U{1,3}{1};

clearvars -except SACCADE INTERVAL SACD Stim U N I T TRIAL FLY GRAND CC Vel

%% Saccade Trial %%
close all ; clc

% velIdx = 5;
% flyIdx = 1;
% trialIdx = 1;
% headIdx = 2;

velIdx = 6;
flyIdx = 12;
trialIdx = 2;
headIdx = 2;

stimVel = Vel(velIdx);


time = TRIAL{flyIdx,velIdx}{trialIdx,headIdx}.Time;
pos  = TRIAL{flyIdx,velIdx}{trialIdx,headIdx}.X(:,1);
vel  = TRIAL{flyIdx,velIdx}{trialIdx,headIdx}.X(:,2);
stim = stimVel*time;


[head.SACD,head.thresh,head.count,head.rate,head.SACDRmv] = SacdDetect(pos,time,2.5,true);
[Saccade,Interval,Stimulus,Error,IntError,matchFlag] = SaccdInter(pos,time,head.SACD, nan, Stim(:,velIdx), true);

FIG = findobj('Type', 'figure');
set(FIG(1:2),'Position',[2 2 6 3]);
sacdColor = [0 0 0.7];
stimColor = [0.5 0.5 0.6];
rawFig = figure ; clf
set(rawFig,'Color','w','Units','inches','Position',[2 2 4 4])
    ax(1) = subplot(2,1,1) ; hold on
    ylabel('Position (°)','FontSize',12)
    h(1) = plot(time,pos,'Color',sacdColor,'LineWidth',1);
    h(2) = plot(time,stim,'Color',stimColor,'LineWidth',1);
    leg = legend(h,'Head','Stimulus','FontSize',12);
    leg.Box = 'off';
ax(2) = subplot(2,1,2) ; hold on
    ylabel('Velocity (°/s)','FontSize',12)
    xlabel('Time (s)','FontSize',12)
    plot(time,vel,'Color',sacdColor,'LineWidth',1)
    plot(time,stimVel*ones(size(time)),'Color',stimColor,'LineWidth',1)
    
set(ax,'FontSize',10,'XLim',[0 4])
linkaxes(ax,'x')

end