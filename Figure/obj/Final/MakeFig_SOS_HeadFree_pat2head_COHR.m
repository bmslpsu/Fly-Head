function [FIG] = MakeFig_SOS_HeadFree_pat2head_COHR()
%% MakeFig_SOS_HeadFree_pat2head_COHR:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     : figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

[Free,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head-free trials', root, 'MultiSelect','on');
Free = cellstr(Free)';

[Fixed,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head-fixed trials', root, 'MultiSelect','on');
Fixed = cellstr(Fixed)';

HeadFree = load(fullfile(root,Free{1}),'GRAND','U','N');
HeadFixed = load(fullfile(root,Fixed{1}),'GRAND','U','N');

figNum = 1;
filename = 'SOS_HeadFree_pat2head_COHR';
xIdx = 1;

FIG = figure (figNum) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 4 3];
FIG.Name = filename;
movegui(FIG,'center')

hold on
ax = gca;
ax.FontSize = 8;
ax.YLim = [0 1];
ax.YLabel.String = 'Coherence';
ax.YLabel.FontSize = 8;
ax.XLim = [0 12];
ax.XLabel.String = 'Frequency (Hz)';
ax.XLabel.FontSize = 8;

PlotPatch(HeadFree.GRAND{5}.Mean{1}{7}(:,xIdx), HeadFree.GRAND{5}.STD{1}{7}(:,xIdx),...
    HeadFree.GRAND{5}.Mean{1}{8}(:,xIdx), 3, HeadFree.N{1,1}, 'b', [0.4 0.4 0.6], 0.5, 2);

[~,h.wing] = PlotPatch(HeadFree.GRAND{8}.Mean{1}{7}(:,xIdx), HeadFree.GRAND{8}.STD{1}{7}(:,xIdx),...
    HeadFree.GRAND{8}.Mean{1}{8}(:,xIdx), 3, HeadFree.N{1,1}, 'r', [0.4 0.4 0.6], 0.5, 2);

[~,h.wing_fixed] = PlotPatch(HeadFixed.GRAND{3}.Mean{1}{7}(:,xIdx), HeadFixed.GRAND{3}.STD{1}{7}(:,xIdx),...
    HeadFixed.GRAND{3}.Mean{1}{8}(:,xIdx), 3, HeadFixed.N{1,1}, 'g', [0.4 0.4 0.6], 0.5, 2);

uistack(h.wing,'top')
uistack(h.wing_fixed,'top')

% errorbar(HeadFree.GRAND{catIdx}.Mean{1}{4}(:,xIdx),HeadFree.GRAND{catIdx}.Mean{1}{9}(:,xIdx),...
%     2*HeadFree.GRAND{catIdx}.STD{1}{9}(:,xIdx),'-or','LineWidth',3)

end