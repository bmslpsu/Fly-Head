function [FIG] = MakeFig_SOS_Walking_HeadFree_Frequency()
%% MakeFig_SOS_HeadFree_pat_head_Time:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'S:\Public\Audrey';

[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select walking trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

Walking = load(fullfile(root,FILE{1}),'TRIAL','GRAND', 'FLY','U','N', 'T');

[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select flying trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';
HeadFree = load(fullfile(root,FILE{1}),'TRIAL','GRAND', 'FLY','U','N', 'T');

%%
figNum = 2;
catIdx = 2;
xIdx = 1;

cpat = [0.4940 0.1840 0.5560];
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


plot(Walking.TRIAL{1}{1,1}.IOFreq, Walking.TRIAL{1}{1,1}.IOMag(:,xIdx),'color', cpat);
PatMean = plot(Walking.TRIAL{1}{1,1}.Fv,Walking.TRIAL{1}{1,1}.Mag(:,xIdx),'color', cpat,'LineWidth',1);

[~, FlyMean] = PlotPatch(HeadFree.GRAND{1,catIdx}.Mean{2}{8}(:,xIdx),HeadFree.GRAND{1,catIdx}.STD{2}{8}(:,xIdx),...
        HeadFree.GRAND{1,catIdx}.Mean{2}{7}(:,xIdx),1,HeadFree.N{1,1},cfly,cfly,0.5,1);
FlyIO = plot(HeadFree.GRAND{1,catIdx}.Mean{2}{10},HeadFree.GRAND{1,catIdx}.Mean{2}{11}(:,xIdx),'color',cfly);

[~, WalkMean] = PlotPatch(Walking.GRAND{1,catIdx}.Mean{2}{8}(:,xIdx),Walking.GRAND{1,catIdx}.STD{2}{8}(:,xIdx),...
        Walking.GRAND{1,catIdx}.Mean{2}{7}(:,xIdx),1,Walking.N{1,1},cwalk,cwalk,0.5,1);
plot(Walking.GRAND{1,catIdx}.Mean{2}{10},Walking.GRAND{1,catIdx}.Mean{2}{11}(:,xIdx),'color',cwalk);

% title('Magtitude vs Frequency of Response')
leg = legend([PatMean WalkMean FlyMean], 'Stimulus', 'Walking Fly', 'Flying Fly');
leg.Box = 'off';


end