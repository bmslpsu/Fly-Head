function [FIG] = MakeFig_SOS_HeadFree_pat2head_BODE()
%% MakeFig_SOS_HeadFree_pat2head_BODE: BODE head position for SOS
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%

root = 'H:\DATA\Rigid_Data\';

% Select files
[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

HeadFree = load(fullfile(PATH,FILE{1}),'TRIAL','GRAND','U','N');

%%
figNum = 1;
catIdx = 5;
xIdx = 1;

filename = 'SOS_HeadFree_pat2head_BODE';

FIG = figure (figNum); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 4 4];
FIG.Name = filename;
movegui(FIG,'center')
hold on

ax1 = subplot(2,1,1) ; hold on
    ax1.FontSize = 8;
    ax1.XLim = [0 10];
    ax1.YLim = [0 0.3];
    ax1.YLim = [0 1];
    ax1.XTickLabel = '';
    ax1.XLabel.FontSize = 8;
    ax1.XLabel.Color = 'none';
 	ax1.YLabel.String = ['Gain (' char(176) '/' char(176) ')'];
    ax1.YLabel.FontSize = ax1.XLabel.FontSize;
    
%     errorbar(HeadFree.GRAND{1,catIdx}.Mean{2}{4},HeadFree.GRAND{1,catIdx}.Mean{2}{5}(:,xIdx),...
%         2*HeadFree.GRAND{1,catIdx}.STD{2}{5}(:,xIdx),'-ob','LineWidth',3);
    
    [~,h.gain] = PlotPatch(HeadFree.GRAND{1,catIdx}.Mean{2}{5}(:,xIdx), HeadFree.GRAND{1,catIdx}.STD{2}{5}(:,xIdx), ...
        HeadFree.GRAND{1,catIdx}.Mean{2}{4}, 1, 1, 'b', [0.4 0.4 0.6], 0.5, 2);
	h.gain.Marker = '.';
    h.gain.MarkerSize = 20;
    
ax2 = subplot(2,1,2) ; hold on
    ax2.FontSize = ax1.FontSize;
    ax2.XLim = ax1.XLim;
    ax2.YLim = [-300 180];
   	ax2.XLabel.String = 'Frequency (Hz)';
    ax2.XLabel.FontSize = 8;
    ax2.XLabel.Color = 'k';
 	ax2.YLabel.String = ['Phase (' char(176) ')'];
    ax2.YLabel.FontSize = ax1.XLabel.FontSize;
    ax2.YTick = -300:60:60;
    
%     errorbar(HeadFree.GRAND{1,catIdx}.Mean{2}{4},rad2deg(HeadFree.GRAND{1,catIdx}.CircMean{7}{6}(:,xIdx)),...
%         2*rad2deg(HeadFree.GRAND{1,catIdx}.CircSTD{7}{6}(:,xIdx)),'-ob','LineWidth',3);
    
    [~,h.phase] = PlotPatch(rad2deg(HeadFree.GRAND{1,catIdx}.CircMean{7}{6}(:,xIdx)), rad2deg(HeadFree.GRAND{1,catIdx}.CircSTD{7}{6}(:,xIdx)), ...
        HeadFree.GRAND{1,catIdx}.Mean{2}{4}, 1, 1, 'b', [0.4 0.4 0.6], 0.5, 2);
    h.phase.Marker = '.';
    h.phase.MarkerSize = 20;
    
    plot(ax1.XLim,[0 0],'--k','LineWidth',1);


end