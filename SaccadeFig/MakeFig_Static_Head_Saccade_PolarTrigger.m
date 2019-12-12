function [FIG] = MakeFig_Static_Head_Saccade_PolarTrigger()
%% MakeFig_Static_Head_Saccade_PolarTrigger:
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
FIG.Position = [2 2 3 3];
FIG.Name = 'Saccade Trigger';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
clear ax

% wave_Idx = SACD.Head.Wave==30;;

ax(1) = subplot(1,1,1,polaraxes); grid off ; axis tight
    h(1) = polarhistogram(deg2rad(SACD.Head.StartPos),100,'FaceColor','g','FaceAlpha',0.9,'Normalization','Probability'); hold on
    h(2) = polarhistogram(deg2rad(SACD.Head.EndPos),100,'FaceColor','r','FaceAlpha',0.9,'Normalization','Probability');;
    ax(1).ThetaAxis.Label.String = ['Head Position (' char(176) ')'];

set(h,'EdgeColor','none')
set(ax,'FontSize',8);
set(ax,'Color','w');
set(ax,'ThetaLim',[-20 20]);
% set(ax,'RLim',[0 55]);
set(ax,'ThetaTick',-20:10:20);
set(ax,'ThetaZeroLocation','top');
leg = legend('Start','End','Location','North');
leg.Location = 'northwest';
leg.Box = 'off';

end