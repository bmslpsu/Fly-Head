function [FIG] = MakeFig_Static_HeadFree_head2wing_COHR()
%% MakeFig_Static_HeadFree_head2wing_COHR:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     : figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

[Static,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head-free trials', root, 'MultiSelect','on');
Static = cellstr(Static)';

[Sos,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head-free trials', root, 'MultiSelect','on');
Sos = cellstr(Sos)';

STATIC = load(fullfile(root,Static{1}),'GRAND','U','N');
SOS = load(fullfile(root,Sos{1}),'GRAND','U','N');

%%
filename = 'Static_HeadFree_head2wing_COHR';
xIdx = 1;
head2wing_stat = 5;
freqIdx = 6;
head2wing_sos = 7;


FIG = figure (1) ; clf
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
ax.XLim = [0 15];
ax.XLabel.String = 'Frequency (Hz)';
ax.XLabel.FontSize = 8;

[~,h.stat] = PlotPatch(STATIC.GRAND{freqIdx,head2wing_stat}.Mean{1}{7}(:,xIdx), STATIC.GRAND{freqIdx,head2wing_stat}.STD{1}{7}(:,xIdx),...
    STATIC.GRAND{freqIdx,head2wing_stat}.Mean{1}{8}(:,xIdx), 1, STATIC.N{1,1}, [0 0.5 0.7], [0.4 0.4 0.6], 0.5, 2);

[~,h.sos] = PlotPatch(SOS.GRAND{1,head2wing_sos}.Mean{1}{7}(:,xIdx), SOS.GRAND{1,head2wing_sos}.STD{1}{7}(:,xIdx),...
    SOS.GRAND{1,head2wing_sos}.Mean{1}{8}(:,xIdx), 1, SOS.N{1,1}, [0.4 0.7 0.2], [0.4 0.4 0.6], 0.5, 2);

uistack([h.sos h.stat], 'top')
leg = legend([h.sos h.stat], 'SOS', 'Static');
leg.Box = 'off';

% errorbar(HeadFree.GRAND{catIdx}.Mean{1}{4}(:,xIdx),HeadFree.GRAND{catIdx}.Mean{1}{9}(:,xIdx),...
%     2*HeadFree.GRAND{catIdx}.STD{1}{9}(:,xIdx),'-or','LineWidth',3)

end