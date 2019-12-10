function [FIG] = MakeFig_Ramp_Head_Saccade_PolarTrigger()
%% MakeFig_Ramp_Head_Saccade_PolarTrigger:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

load(fullfile(root,FILE{1}),'SACD','U','N');

clearvars -except SACCADE INTERVAL SACD Stim U N I TRIAL FLY GRAND CC Vel

%% Saccade Trigger Polar Plot %%
FIG = figure (1) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 6 3];
FIG.Name = 'Saccade Trigger';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
clear ax

match = -1;
mIdx = SACD.Head.Match==match;
MATCH = SACD.Head(mIdx,:);

ax(1) = subplot(1,2,1,polaraxes); grid off ; axis tight
    h(1) = polarhistogram(deg2rad(MATCH.StartPos(MATCH.Dir==1)),100,'FaceColor','g','FaceAlpha',0.9,'Normalization','Probability'); hold on
    h(2) = polarhistogram(deg2rad(MATCH.EndPos(MATCH.Dir==1)),100,'FaceColor','r','FaceAlpha',0.9,'Normalization','Probability');
    title('CW Stimulus Direction   ===>')
    ax(1).ThetaAxis.Label.String = ['Head Position (' char(176) ')'];
ax(2) = subplot(1,2,2,polaraxes); grid off ; axis tight
    h(3) = polarhistogram(deg2rad(MATCH.StartPos(MATCH.Dir==-1)),100,'FaceColor','g','FaceAlpha',0.9,'Normalization','Probability'); hold on
    h(4) = polarhistogram(deg2rad(MATCH.EndPos(MATCH.Dir==-1)),100,'FaceColor','r','FaceAlpha',0.9,'Normalization','Probability');
    title('CCW Stimulus Direction  <===')
    ax(2).ThetaAxis.Label.String = ['Head Position (' char(176) ')'];

set(h,'EdgeColor','none')
set(ax,'FontSize',8);
set(ax,'Color','w');
set(ax,'ThetaLim',[-20 20]);
% set(ax,'RLim',[0 300]);
set(ax,'ThetaDir','clockwise')
set(ax,'ThetaTick',-20:10:20);
set(ax,'ThetaZeroLocation','top');
leg = legend('Start','End','Location','North');
leg.Location = 'northwest';
leg.Box = 'off';

end