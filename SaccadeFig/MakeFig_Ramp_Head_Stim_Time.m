function [FIG] = MakeFig_Ramp_Head_Stim_Time()
%% MakeFig_Ramp_Head_Stim_Time:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';
[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

% load(fullfile(root,FILE{1}),'SACCADE','INTERVAL','SACD','Stim','U','N','I','TRIAL','FLY','GRAND');
load(fullfile(root,FILE{1}),'SACD','Stim','U','N','TRIAL','FLY','GRAND');

clearvars -except SACCADE INTERVAL SACD Stim U N I TRIAL FLY GRAND

CC = repmat({[1 0 0],[0 1 0],[0 0 1]},1,2);
Vel = 3.75*U{1,3}{1};

%% Saccade Removed Time Domain %%
FIG = figure (1) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 4 2];
FIG.Name = 'Removed Saccades';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
hold on
clear h ax

catIdx = [1 2];
xIdx = 2;
velIdx = 5;
flyIdx = 1;
trialIdx = 1;

ax(1) = subplot(1,1,1) ; hold on
ax(1).Title.String = [num2str(Vel(velIdx)) ' (' char(176) '/s)'];
ax(1).XLabel.String = 'Time(s)';
ax(1).YLabel.String = ['Velocity (' char(176) '/s)'];

h.stim = plot(TRIAL{flyIdx,velIdx}{trialIdx,catIdx(1)}.Time, TRIAL{flyIdx,velIdx}{trialIdx,catIdx(1)}.X(:,xIdx), ...
    'Color', 'k', 'LineWidth',1);

h.head = plot(TRIAL{flyIdx,velIdx}{trialIdx,catIdx(2)}.Time, TRIAL{flyIdx,velIdx}{trialIdx,catIdx(2)}.X(:,xIdx), ...
    'Color', 0.7*CC{velIdx}, 'LineWidth',1);

leg = legend([h.head h.stim],'Head','Stimulus');
leg.Box = 'off';

set(ax,'FontSize',8)
set(ax,'XLim',[0 10],'YLim',1000*[-1 1])
linkaxes(ax,'x')

end