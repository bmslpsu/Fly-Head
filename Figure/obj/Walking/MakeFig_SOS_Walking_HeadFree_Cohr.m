function [FIG] = MakeFig_SOS_Walking_HeadFree_Cohr()
%% MakeFig_SOS_HeadFree_pat_head_Time:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'E:\Walking_Experiments';

[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select walking trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

Walking = load(fullfile(root,FILE{1}),'TRIAL','GRAND', 'FLY','U','N', 'T');

root = 'S:\Public\Audrey';
[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select flying trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';
HeadFree = load(fullfile(root,FILE{1}),'TRIAL','GRAND', 'FLY','U','N', 'T');

%%
figNum = 3;
catIdx = 4;
xIdx = 1;

cwalk = [0 0.4470 0.7410];
cfly = [0.8500 0.3250 0.0980];

filename = 'SOS_Walking_all_head_Time';

FIG = figure (figNum); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 3 3];
movegui(FIG,'center')
FIG.Name = filename;

ax.L = gca; hold on
ax.L.YColor = 'k';
ax.L.FontSize = 8;
ax.L.YLabel.String = ['Magnitude'];
ax.L.YLabel.FontSize = 10;
% ax.L.YLim =[0 10];
% ax.L.YTick = 15*[-1 0 1];
ax.L.XLabel.String = 'Frequency (Hz)';
ax.L.XLabel.Color = 'k';
ax.L.XLabel.FontSize = 10;
ax.L.XLim = [.5 10];
ax.L.XTick = [1, 3.1, 5.3, 7.4, 9.6];


[~, WalkCohr] = PlotPatch(Walking.GRAND{1,catIdx}.Mean{2}{7}(:,xIdx),Walking.GRAND{1,catIdx}.STD{2}{7}(:,xIdx), ...
    Walking.GRAND{1,catIdx}.Mean{2}{8}(:,xIdx),1,Walking.N{1,1},cwalk,cwalk,0.5,1);
[~, WalkCohrConn] = PlotPatch(Walking.GRAND{1,catIdx}.Mean{2}{9}, Walking.GRAND{1,catIdx}.STD{2}{9},...
    Walking.GRAND{1,catIdx}.Mean{2}{4},1,Walking.N{1,1},cwalk,cwalk,0.5,1);


catIdx = 5;
[~, FlyCohr] = PlotPatch(HeadFree.GRAND{1,catIdx}.Mean{2}{7}(:,xIdx),HeadFree.GRAND{1,catIdx}.STD{2}{7}(:,xIdx), ...
    HeadFree.GRAND{1,catIdx}.Mean{2}{8}(:,xIdx),1,HeadFree.N{1,1},cfly,cfly,0.5,1);
[~,FlyCohrConn] = PlotPatch(HeadFree.GRAND{1,catIdx}.Mean{2}{9}(:,xIdx),HeadFree.GRAND{1,catIdx}.STD{2}{9}(:,xIdx), ...
    HeadFree.GRAND{1,catIdx}.Mean{2}{4},1,HeadFree.N{1,1},cfly,cfly,0.5,1);


%title('Coherence of Fly Response')
leg = legend([WalkCohr FlyCohr], 'Walking Fly', 'Flying Fly')
leg.Box = 'off'

end