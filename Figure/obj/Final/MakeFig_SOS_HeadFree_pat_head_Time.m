function [FIG] = MakeFig_SOS_HeadFree_pat_head_Time()
%% MakeFig_SOS_HeadFree_pat_head_Time:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

HeadFree = load(fullfile(root,FILE{1}),'TRIAL','GRAND','U','N');

figNum = 1;
catIdx = 2;
xIdx = 1;

filename = 'SOS_HeadFree_pat_head_Time';

FIG = figure (figNum); clf
FIG.Color = 'w';
FIG.Position = [100 100 1200 400];
movegui(FIG,'center')
FIG.Name = filename;

yyaxis left
ax.L = gca; hold on
ax.L.YColor = 'k';
ax.L.FontSize = 12;
ax.L.YLabel.String = ['Head (' char(176) ')'];
ax.L.YLabel.FontSize = 14;
ax.L.YLim = 20*[-1 1];
ax.L.XLabel.String = 'Time (s)';
ax.L.XLabel.Color = 'k';
ax.L.XLabel.FontSize = ax.L.YLabel.FontSize;
ax.L.XLim = [0 20];

PlotPatch(HeadFree.GRAND{1,catIdx}.Mean{2}{6}(:,xIdx),HeadFree.GRAND{1,catIdx}.STD{2}{6}(:,xIdx),...
    HeadFree.GRAND{1,catIdx}.Mean{2}{5},3, HeadFree.N{1,1}, 'b', [0.4 0.4 0.6], 0.5, 2);

yyaxis right
ax.R = gca; hold on
ax.R.YColor = 'g';
ax.R.FontSize = ax.L.FontSize;
ax.R.YLabel.String = ['Stimulus (' char(176) ')'];
ax.R.YLabel.FontSize = ax.L.YLabel.FontSize;
ax.R.YLim = 20*[-1 1];
ax.R.XLabel.FontSize = ax.R.YLabel.FontSize;
ax.R.XLim = ax.L.XLim;

plot(HeadFree.TRIAL{1}{2,1}.Time,HeadFree.TRIAL{1}{2,1}.X(:,xIdx),'g','LineWidth',2)

end