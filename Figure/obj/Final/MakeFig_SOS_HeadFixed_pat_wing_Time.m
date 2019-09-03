function [FIG] = MakeFig_SOS_HeadFixed_pat_wing_Time()
%% MakeFig_SOS_HeadFixed_pat_wing_Time:
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
%%
figNum = 1;
catIdx = 2;
xIdx = 1;

filename = 'SOS_HeadFree_pat_wing_Time';

FIG = figure (figNum); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 8 2];
movegui(FIG,'center')
FIG.Name = filename;

ax.R = gca; hold on
ax.R.YColor = 'g';
ax.R.YAxisLocation = 'right';
ax.R.FontSize = 8;
ax.R.YLabel.String = ['(' char(176) ')'];
ax.R.YLabel.FontSize = 8;
ax.R.YLim = 15*[-1 1];
ax.R.YTick = 15*[-1 0 1];
ax.R.XLabel.String = 'Time (s)';
ax.R.XLabel.Color = 'k';
ax.R.XLabel.FontSize = 8;
ax.R.XLim = [0 20];

plot(HeadFree.TRIAL{1}{2,1}.Time,HeadFree.TRIAL{1}{2,1}.X(:,xIdx),'g','LineWidth',2)

ax.L = axes; hold on
ax.L.YColor = 'k';
ax.L.FontSize = 8;
ax.L.YLabel.String = '\Delta WBA (V)';
ax.L.YLabel.FontSize = 8;
ax.L.YLim = 2*[-1 1];
ax.L.XLabel.FontSize = 8;
ax.L.XLim = ax.R.XLim;
ax.L.Color = 'none';
ax.L.Position = ax.R.Position;
ax.L.XColor = 'none';

PlotPatch(HeadFree.GRAND{1,catIdx}.Mean{2}{6}(:,xIdx),HeadFree.GRAND{1,catIdx}.STD{2}{6}(:,xIdx),...
    HeadFree.GRAND{1,catIdx}.Mean{2}{5},2, HeadFree.N{1,1}, 'r', [0.4 0.4 0.6], 0.5, 2);

end