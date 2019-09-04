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

load(fullfile(root,FILE{1}),'SACCADE','INTERVAL','SACD','Stim','U','N','I','TRIAL','FLY','GRAND');

CC = repmat({[1 0 0],[0 1 0],[0 0 1]},1,2);

Vel = 3.75*U{1,3}{1};

clearvars -except SACCADE INTERVAL SACD Stim U N I TRIAL FLY GRAND CC Vel

%% Saccade Trial %%
FIG = figure (1) ; clf
FIG.Units = 'inches';
FIG.Position = 0.75*[2 2 3 2];
FIG.Name = 'Saccade Position';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
hold on

velIdx = 1;
flyIdx = 1;
trialIdx = 1;
headIdx = 2;

time = TRIAL{flyIdx,velIdx}{trialIdx,headIdx}.Time;
pos  = TRIAL{flyIdx,velIdx}{trialIdx,headIdx}.X(:,1);
% pos = pos - pos(1);

plot(time,Stim(:,velIdx),'-','Color',[0.5 0.5 0.5],'LineWidth',1)
plot(time,pos,'-','Color',[0.5 0.1 0.8],'LineWidth',1.5)
xlim([0 3])
xlabel('Time (s)')
ylabel(['Angle (' char(176) ')'])
leg = legend('Stimulus','Head','location','northwest');
leg.Box = 'off';

end