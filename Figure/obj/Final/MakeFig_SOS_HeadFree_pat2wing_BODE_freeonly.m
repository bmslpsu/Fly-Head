function [FIG] = MakeFig_SOS_HeadFree_pat2wing_BODE_freeonly()
%% MakeFig_SOS_HeadFree_pat2wing_BODE_freeonly:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%

root = 'H:\DATA\Rigid_Data\';

% Select files
[Free,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head-free data', root, 'MultiSelect','off');
Free = cellstr(Free)';

HeadFree = load(fullfile(root,Free{1}),'GRAND','U','N');

%%
figNum = 1;
xIdx = 1;

filename = 'SOS_HeadFree_pat2wing_BODE';

stim2wing = 8;
err2wing = 6;

FIG = figure (figNum); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 4 4];
FIG.Name = filename;
movegui(FIG,'center')
hold on

Phase_2 = HeadFree.GRAND{1,err2wing}.CircMean{7}{6}(:,xIdx);
cnd = (Phase_2>2) & ((1:length(Phase_2))>1)';
Phase_2(cnd) = Phase_2(cnd) - 2*pi;

ax1 = subplot(2,1,1) ; hold on
    ax1.FontSize = 8;
    ax1.XLim = [0 10];
	ax1.YLim = [0 0.35];
    ax1.XTickLabel = '';
    ax1.XLabel.FontSize = 8;
    ax1.XLabel.Color = 'none';
 	ax1.YLabel.String = ['Gain (' char(176) '/' char(176) ')'];
    ax1.YLabel.FontSize = 8;
    
    [~,h.gain_1] = PlotPatch(HeadFree.GRAND{1,stim2wing}.Mean{2}{5}(:,xIdx), HeadFree.GRAND{1,stim2wing}.STD{2}{5}(:,xIdx), ...
        HeadFree.GRAND{1,stim2wing}.Mean{2}{4}, 1, 1, 'r', [0.4 0.4 0.6], 0.5, 2);
	h.gain_1.Marker = '.';
    h.gain_1.MarkerSize = 20;
    
  	[~,h.gain_2] = PlotPatch(HeadFree.GRAND{1,err2wing}.Mean{2}{5}(:,xIdx), HeadFree.GRAND{1,err2wing}.STD{2}{5}(:,xIdx), ...
        HeadFree.GRAND{1,err2wing}.Mean{2}{4}, 1, 1, 'c', [0.4 0.4 0.6], 0.5, 2);
	h.gain_2.Marker = '.';
    h.gain_2.MarkerSize = 20;
    
    uistack(h.gain_2,'top')
    uistack(h.gain_1,'top')
	
ax2 = subplot(2,1,2) ; hold on
    ax2.FontSize = ax1.FontSize;
    ax2.XLim = ax1.XLim;
    ax2.YLim = [-300 150];
   	ax2.XLabel.String = 'Frequency (Hz)';
    ax2.XLabel.FontSize = 8;
    ax2.XLabel.Color = 'k';
 	ax2.YLabel.String = ['Phase (' char(176) ')'];
    ax2.YLabel.FontSize = 8;
    ax2.YTick = -300:60:120;
    
    [~,h.phase_1] = PlotPatch(rad2deg(HeadFree.GRAND{1,stim2wing}.CircMean{7}{6}(:,xIdx)), rad2deg(HeadFree.GRAND{1,stim2wing}.CircSTD{7}{6}(:,xIdx)), ...
        HeadFree.GRAND{1,stim2wing}.Mean{2}{4}, 1, 1, 'r', [0.4 0.4 0.6], 0.5, 2);
    h.phase_1.Marker = '.';
    h.phase_1.MarkerSize = 20;
    
 	[~,h.phase_2] = PlotPatch(rad2deg(Phase_2), rad2deg(HeadFree.GRAND{1,err2wing}.CircSTD{7}{6}(:,xIdx)), ...
        HeadFree.GRAND{1,err2wing}.Mean{2}{4}, 1, 1, 'c', [0.4 0.4 0.6], 0.5, 2);
    h.phase_2.Marker = '.';
    h.phase_2.MarkerSize = 20;
    
    uistack(h.phase_2,'top')
    uistack(h.phase_1,'top')
    
    plot(ax1.XLim,[0 0],'--k','LineWidth',1);
    
    leg = legend([h.phase_1,h.phase_2],'Free: Ref2Wing','Free: Err2Wing');
    leg.Box = 'off';


end