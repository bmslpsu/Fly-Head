function [FIG] = MakeFig_SOS_HeadFree_pat2head_BODE()
%% MakeFig_SOS_HeadFree_pat2head_BODE: BODE head position for SOS
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

% Select files
[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

HeadFree = load(fullfile(root,FILE{1}),'TRIAL','FLY','GRAND','U','N');

figNum = 1;
catIdx = 5;
xIdx = 1;

filename = 'SOS_HeadFree_pat2head_BODE';

FIG = figure (figNum); clf
FIG.Color = 'w';
FIG.Position = [100 100 600 600];
FIG.Name = filename;
movegui(FIG,'center')
hold on


ax1 = subplot(2,1,1) ; hold on
    ax1.FontSize = 12;
    ax1.XLim = [0 10]; 
    ax1.YLim = [0 1];
    ax1.XTickLabel = '';
% 	ax1.XLabel.String = 'Frequency (Hz)';
    ax1.XLabel.FontSize = 14;
    ax1.XLabel.Color = 'w';
 	ax1.YLabel.String = ['Gain (' char(176) '/' char(176) ')'];
    ax1.YLabel.FontSize = ax1.XLabel.FontSize;
    
    errorbar(HeadFree.GRAND{1,catIdx}.Mean{2}{4},HeadFree.GRAND{1,catIdx}.Mean{2}{5}(:,xIdx),...
        2*HeadFree.GRAND{1,catIdx}.STD{2}{5}(:,xIdx),'-ob','LineWidth',3);
    
ax2 = subplot(2,1,2) ; hold on
    ax2.FontSize = ax1.FontSize;
    ax2.XLim = ax1.XLim;
    ax2.YLim = [-300 60];
   	ax2.XLabel.String = 'Frequency (Hz)';
    ax2.XLabel.FontSize = 14;
    ax2.XLabel.Color = 'k';
 	ax2.YLabel.String = ['Phase (' char(176) ')'];
    ax2.YLabel.FontSize = ax1.XLabel.FontSize;
    ax2.YTick = -300:60:60;
    
    errorbar(HeadFree.GRAND{1,catIdx}.Mean{2}{4},rad2deg(HeadFree.GRAND{1,catIdx}.CircMean{9}{6}(:,xIdx)),...
        rad2deg(2*HeadFree.GRAND{1,catIdx}.CircSTD{9}{6}(:,xIdx)),'-ob','LineWidth',3);
    
    plot(ax1.XLim,[0 0],'--g','LineWidth',2);


end