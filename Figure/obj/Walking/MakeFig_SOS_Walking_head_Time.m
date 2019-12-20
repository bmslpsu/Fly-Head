function [FIG] = MakeFig_SOS_Walking_head_Time()
%% MakeFig_SOS_HeadFree_pat_head_Time:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'S:\Public\Audrey\Walking SOS';

[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select Walking Data', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

Walking = load(fullfile(root,FILE{1}),'TRIAL','GRAND', 'FLY','U','N', 'T');

[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select Flight Data', root, 'MultiSelect','off');
FILE = cellstr(FILE)';
HeadFree = load(fullfile(root,FILE{1}),'TRIAL','GRAND', 'FLY','U','N', 'T');
%%
figNum = 1;
catIdx = 2;
xIdx = 1;

cpat = [0.4940 0.1840 0.5560];
cwalk = [0 0.4470 0.7410];
cfly = [0.8500 0.3250 0.0980];

filename = 'SOS_Walking_all_head_Time';

FIG = figure (figNum); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 6 3];
movegui(FIG,'center')
FIG.Name = filename;

ax.L = gca; hold on
ax.L.YColor = 'k';
ax.L.FontSize = 8;
ax.L.YLabel.String = ['Head (' char(176) ')'];
ax.L.YLabel.FontSize = 10;
ax.L.YLim = 20*[-1 1];
ax.L.YTick = 20*[-1 0 1];
ax.L.XLabel.String = 'Time (s)';
ax.L.XLabel.Color = 'k';
ax.L.XLabel.FontSize = 10;
ax.L.XLim = [0 20];


Pat = plot(Walking.TRIAL{1}{1,1}.Time,Walking.TRIAL{1}{1,1}.X(:,xIdx),'color', cpat,'LineWidth',1);
 
[~, FlyMag] = PlotPatch(HeadFree.GRAND{1,catIdx}.Mean{2}{6}(:,xIdx),HeadFree.GRAND{1,catIdx}.STD{2}{6}(:,xIdx),...
        HeadFree.GRAND{1,catIdx}.Mean{2}{5},1, HeadFree.N{1,1},cfly,cfly,0.3,1);
    
[~, WalkMag] = PlotPatch(Walking.GRAND{1,catIdx}.Mean{2}{6}(:,xIdx),Walking.GRAND{1,catIdx}.STD{2}{6}(:,xIdx),...
        Walking.GRAND{1,catIdx}.Mean{2}{5},1,Walking.N{1,1},cwalk,cwalk,0.3,1);
    

leg=legend([Pat WalkMag FlyMag], 'Pattern', 'Walking Fly', 'Flying Fly')
leg.Box = 'off'

 

end