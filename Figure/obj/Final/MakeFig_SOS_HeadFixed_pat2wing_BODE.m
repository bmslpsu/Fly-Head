function [FIG] = MakeFig_SOS_HeadFixed_pat2wing_BODE()
%% MakeFig_SOS_HeadFixed_pat2wing_BODE:
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

HeadFixed = load(fullfile(root,FILE{1}),'TRIAL','FLY','GRAND','U','N');

catIdx = 3;
xIdx = 1;

filename = 'SOS_HeadFixed_pat2wing_BODE';
%%
FIG = figure (1); clf
FIG.Color = 'w';
FIG.Position = [100 100 600 600];
FIG.Name = filename;
movegui(FIG,'center')
hold on

ax1 = subplot(2,1,1) ; hold on
    ax1.FontSize = 12;
    ax1.XLim = [0 10]; 
    ax1.YLim = [0 0.3];
    ax1.XTickLabel = '';
% 	ax1.XLabel.String = 'Frequency (Hz)';
    ax1.XLabel.FontSize = 14;
    ax1.XLabel.Color = 'w';
 	ax1.YLabel.String = ['Gain (V/' char(176) ')'];
    ax1.YLabel.FontSize = ax1.XLabel.FontSize;
    
    errorbar(HeadFixed.GRAND{1,catIdx}.Mean{2}{4},HeadFixed.GRAND{1,catIdx}.Mean{2}{5}(:,xIdx),...
        2*HeadFixed.GRAND{1,catIdx}.STD{2}{5}(:,xIdx),'-oc','LineWidth',3);
    
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
    
    GAIN = rad2deg(HeadFixed.GRAND{1,catIdx}.CircMean{9}{6}(:,xIdx));
    GAIN(GAIN>60) = GAIN(GAIN>60) - 360;
    
    errorbar(HeadFixed.GRAND{1,catIdx}.Mean{2}{4},GAIN,...
        rad2deg(2*HeadFixed.GRAND{1,catIdx}.CircSTD{9}{6}(:,xIdx)),'-oc','LineWidth',3);
    
    plot(ax1.XLim,[0 0],'--g','LineWidth',2);
%%
FIG = figure (2); clf
FIG.Color = 'w';
FIG.Position = [100 100 600 600];
FIG.Name = filename;
movegui(FIG,'center')
hold on

ax1 = subplot(2,1,1) ; hold on
    ax1.FontSize = 12;
    ax1.XLim = [0 10]; 
    ax1.YLim = [0 0.3];
    ax1.XTickLabel = '';
% 	ax1.XLabel.String = 'Frequency (Hz)';
    ax1.XLabel.FontSize = 14;
    ax1.XLabel.Color = 'w';
 	ax1.YLabel.String = ['Gain (V/' char(176) ')'];
    ax1.YLabel.FontSize = ax1.XLabel.FontSize;
    
    errorbar(HeadFixed.GRAND{1,catIdx}.Mean{2}{4},HeadFixed.GRAND{1,catIdx}.Mean{2}{5}(:,xIdx),...
        2*HeadFixed.GRAND{1,catIdx}.STD{2}{5}(:,xIdx),'-oc','LineWidth',3);
    
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
    
    GAIN = rad2deg(HeadFixed.GRAND{1,catIdx}.CircMean{9}{6}(:,xIdx));
    GAIN(GAIN>60) = GAIN(GAIN>60) - 360;
    
    errorbar(HeadFixed.GRAND{1,catIdx}.Mean{2}{4},GAIN,...
        rad2deg(2*HeadFixed.GRAND{1,catIdx}.CircSTD{9}{6}(:,xIdx)),'-oc','LineWidth',3);
    
    plot(ax1.XLim,[0 0],'--g','LineWidth',2);
    
end