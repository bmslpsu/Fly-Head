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

velIdx = 5;
flyIdx = 1;
trialIdx = 1;
headIdx = 2;

time = TRIAL{flyIdx,velIdx}{trialIdx,headIdx}.Time;
pos  = TRIAL{flyIdx,velIdx}{trialIdx,headIdx}.X(:,1);

[head.SACD,head.thresh,head.count,head.rate,head.SACDRmv] = SacdDetect(pos,time,2.5,true);
[Saccade,Interval,Stimulus,Error,IntError,matchFlag] = SaccdInter(pos,time,head.SACD, nan, Stim(:,velIdx), true);

FIG = findobj('Type', 'figure');
set(FIG(1:2),'Position',[2 2 6 3]);

end