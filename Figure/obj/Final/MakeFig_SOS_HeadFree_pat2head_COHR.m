function [FIG] = MakeFig_SOS_HeadFree_pat2head_COHR()
%% MakeFig_SOS_HeadFree_pat2head_COHR:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     : figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

[FILES,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','on');
FILES = cellstr(FILES)';

HeadFree = load(fullfile(root,FILES{1}),'GRAND','U','N');

figNum = 1;
filename = 'SOS_HeadFree_pat2head_COHR';
catIdx = 3;
xIdx = 1;

FIG = figure (figNum) ; clf
FIG.Color = 'w';
FIG.Position = [100 100 680 500];
FIG.Name = filename;
movegui(FIG,'center')

hold on
ax = gca;
ax.FontSize = 12;
ax.YLim = [0 1];
ax.YLabel.String = 'Coherence';
ax.YLabel.FontSize = 14;
ax.XLim = [0 12];
ax.XLabel.String = 'Frequency (Hz)';
ax.XLabel.FontSize = ax.YLabel.FontSize;

PlotPatch(HeadFree.GRAND{catIdx}.Mean{1}{7}(:,xIdx), HeadFree.GRAND{catIdx}.STD{1}{7}(:,xIdx),...
    HeadFree.GRAND{catIdx}.Mean{1}{8}(:,xIdx), 3, HeadFree.N{1,1}, 'k', [0.4 0.4 0.6], 0.5, 2);

% errorbar(HeadFree.GRAND{catIdx}.Mean{1}{4}(:,xIdx),HeadFree.GRAND{catIdx}.Mean{1}{9}(:,xIdx),...
%     2*HeadFree.GRAND{catIdx}.STD{1}{9}(:,xIdx),'-or','LineWidth',3)

end