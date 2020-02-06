function [FIG] = MakeFig_SOS_HeadFree_head2wing_BODE()
%% MakeFig_SOS_HeadFree_head2wing_BODE:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%
root = 'H:\DATA\Rigid_Data\';

% Select files
[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

HeadFree = load(fullfile(root,FILE{1}),'TRIAL','FLY','GRAND','U','N');
%%
figNum = 1;
catIdx = 7;
xIdx = 1;
filename = 'SOS_HeadFree_head2wing_BODE';

IOFreq = HeadFree.GRAND{1,catIdx}.Mean{2}{4};
Gain = HeadFree.GRAND{1,catIdx}.Mean{2}{5}(:,xIdx);
Gain_STD = HeadFree.GRAND{1,catIdx}.STD{2}{5}(:,xIdx);
Phase = rad2deg(HeadFree.GRAND{1,catIdx}.CircMean{7}{6}(:,xIdx));
Phase(Phase>60) = Phase(Phase>60) - 360;
Phase_STD = rad2deg(HeadFree.GRAND{1,catIdx}.CircSTD{7}{6}(:,xIdx));
TimeDelay = 1000*(Phase./360).*(1./IOFreq);
TimeDelay_STD = 1000*(Phase_STD./360).*(1./IOFreq);


FIG = figure (figNum); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 4 4];
FIG.Name = filename;
movegui(FIG,'center')
hold on

ax1 = subplot(3,1,1) ; hold on
    ax1.FontSize = 8;
    ax1.XLim = [0 10]; 
    ax1.YLim = [0 0.3];
    ax1.XTickLabel = '';
    ax1.XLabel.FontSize = 8;
    ax1.XLabel.Color = 'w';
 	ax1.YLabel.String = ['Gain (V/' char(176) ')'];
    ax1.YLabel.FontSize = ax1.XLabel.FontSize;
    ax1.XTickLabels = '';
        
    [~,h.gain] = PlotPatch(Gain , Gain_STD, IOFreq, 1, HeadFree.N{1,1}, 'k', [0.4 0.4 0.6], 0.5, 2);
    h.gain.Marker = '.';
    h.gain.MarkerSize = 15;
    
ax2 = subplot(3,1,2) ; hold on
    ax2.FontSize = 8;
    ax2.XLim = ax1.XLim;
    ax2.YLim = [-180 0];
    ax2.XLabel.FontSize = 8;
    ax2.XLabel.Color = 'k';
 	ax2.YLabel.String = ['Phase (' char(176) ')'];
    ax2.YLabel.FontSize = ax1.XLabel.FontSize;
    ax2.YTick = -180:60:0;
    ax2.XTickLabels = '';
    
	[~,h.phase] = PlotPatch(Phase , Phase_STD, IOFreq, 1, HeadFree.N{1,1}, 'k', [0.4 0.4 0.6], 0.5, 2);
    h.phase.Marker = '.';
    h.phase.MarkerSize = 15;
        
ax3 = subplot(3,1,3) ; hold on
    ax3.FontSize = 8;
    ax3.XLim = ax1.XLim;
   	ax3.XLabel.String = 'Frequency (Hz)';
    ax3.XLabel.FontSize = 8;
    ax3.XLabel.Color = 'k';
 	ax3.YLabel.String = 'Time Difference (ms)';
    ax3.YLabel.FontSize = ax1.XLabel.FontSize;
    ax3.YLim = [-80 0];
    ax3.YTick = -80:20:0;
    
	[~,h.time] = PlotPatch(TimeDelay , TimeDelay_STD, IOFreq, 1, HeadFree.N{1,1}, 'k', [0.4 0.4 0.6], 0.5, 2);
    h.time.Marker = '.';
    h.time.MarkerSize = 15;
    end