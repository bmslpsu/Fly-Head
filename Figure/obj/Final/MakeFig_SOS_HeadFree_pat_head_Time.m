function [FIG] = MakeFig_SOS_HeadFree_pat_head_Time()
%% MakeFig_SOS_HeadFree_pat_head_Time:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%

root = 'H:\DATA\Rigid_Data\';

[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

HeadFree = load(fullfile(root,FILE{1}),'TRIAL','GRAND','U','N');

%%
figNum = 1;
catIdx = 2;
xIdx = 1;

filename = 'SOS_HeadFree_pat_head_Time';

FIG = figure (figNum); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 8 2];
movegui(FIG,'center')
FIG.Name = filename;

ax.L = gca; hold on
ax.L.YColor = 'k';
ax.L.FontSize = 8;
ax.L.YLabel.String = ['Head (' char(176) ')'];
ax.L.YLabel.FontSize = 8;
ax.L.YLim = 15*[-1 1];
ax.L.YTick = 15*[-1 0 1];
ax.L.XLabel.String = 'Time (s)';
ax.L.XLabel.Color = 'k';
ax.L.XLabel.FontSize = 8;
ax.L.XLim = [0 20];

plot(HeadFree.TRIAL{1}{2,1}.Time,HeadFree.TRIAL{1}{2,1}.X(:,xIdx),'k','LineWidth',0.5)

PlotPatch(HeadFree.GRAND{1,catIdx}.Mean{1}{6}(:,xIdx),HeadFree.GRAND{1,catIdx}.STD{2}{6}(:,xIdx),...
    HeadFree.GRAND{1,catIdx}.Mean{2}{5}, 1, HeadFree.N{1,1}, 'b', [0.4 0.4 0.6], 0.5, 1);

end