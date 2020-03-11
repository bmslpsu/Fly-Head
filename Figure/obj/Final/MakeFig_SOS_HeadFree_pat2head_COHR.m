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
%%
filename = 'SOS_HeadFree_pat2head_COHR';
xIdx = 1;
pat2head = 5;
pat2wing = 8;
fixedwing = 3;
head2wing = 7;

FIG = figure (1) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 3 3];
FIG.Name = filename;
movegui(FIG,'center')

hold on
ax = gca;
ax.FontSize = 8;
ax.YLim = [0 1];
ax.YLabel.String = 'Coherence';
ax.YLabel.FontSize = 8;
ax.XLim = [0 10];
ax.XLabel.String = 'Frequency (Hz)';
ax.XLabel.FontSize = 8;

PlotPatch(HeadFree.GRAND{pat2head}.Mean{1}{7}(:,xIdx), HeadFree.GRAND{pat2head}.STD{1}{7}(:,xIdx),...
    HeadFree.GRAND{pat2head}.Mean{1}{8}(:,xIdx), 1, HeadFree.N{1,1}, 'b', [0.4 0.4 0.6], 0.5, 2);

[~,h.wing] = PlotPatch(HeadFree.GRAND{pat2wing}.Mean{1}{7}(:,xIdx), HeadFree.GRAND{pat2wing}.STD{1}{7}(:,xIdx),...
    HeadFree.GRAND{pat2wing}.Mean{1}{8}(:,xIdx), 1, HeadFree.N{1,1}, 'r', [0.4 0.4 0.6], 0.5, 2);

[~,h.wing_fixed] = PlotPatch(HeadFixed.GRAND{fixedwing}.Mean{1}{7}(:,xIdx), HeadFixed.GRAND{fixedwing}.STD{1}{7}(:,xIdx),...
    HeadFixed.GRAND{fixedwing}.Mean{1}{8}(:,xIdx), 1, HeadFixed.N{1,1}, [0.4 0 0.8], [0.4 0.4 0.6], 0.5, 2);

[~,h.headwing] = PlotPatch(HeadFree.GRAND{head2wing}.Mean{1}{7}(:,xIdx), HeadFree.GRAND{head2wing}.STD{1}{7}(:,xIdx),...
    HeadFree.GRAND{head2wing}.Mean{1}{8}(:,xIdx), 1, HeadFree.N{1,1}, [0 0.8 0.2], [0.4 0.4 0.6], 0.5, 2);

uistack(h.wing,'top')
uistack(h.wing_fixed,'top')

errorbar(HeadFree.GRAND{pat2head}.Mean{1}{4}(:,xIdx),HeadFree.GRAND{pat2head}.Mean{1}{9}(:,xIdx),...
    1*HeadFree.GRAND{pat2head}.STD{1}{9}(:,xIdx)./sqrt(HeadFree.N{1,1}),'-b','LineWidth',2)

errorbar(HeadFree.GRAND{pat2wing}.Mean{1}{4}(:,xIdx),HeadFree.GRAND{pat2wing}.Mean{1}{9}(:,xIdx),...
    1*HeadFree.GRAND{pat2wing}.STD{1}{9}(:,xIdx)./sqrt(HeadFree.N{1,1}),'-r','LineWidth',2)

errorbar(HeadFixed.GRAND{fixedwing}.Mean{1}{4}(:,xIdx),HeadFixed.GRAND{fixedwing}.Mean{1}{9}(:,xIdx),...
    1*HeadFixed.GRAND{fixedwing}.STD{1}{9}(:,xIdx)./sqrt(HeadFixed.N{1,1}),'Color',[0.4 0 0.8],'LineWidth',2)

errorbar(HeadFree.GRAND{head2wing}.Mean{1}{4}(:,xIdx),HeadFree.GRAND{head2wing}.Mean{1}{9}(:,xIdx),...
    1*HeadFree.GRAND{head2wing}.STD{1}{9}(:,xIdx)./sqrt(HeadFree.N{1,1}), 'Color', [0 0.8 0.2], 'LineWidth',2)

set(ax,'LineWidth',1.5,'FontSize',8)
end