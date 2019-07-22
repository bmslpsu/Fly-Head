function [FIG] = MakeFig_SOS_HeadFree_head_wing_TimeDiff()
%% MakeFig_SOS_HeadFree_head_wing_TimeDiff: 
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

HeadFree = load(fullfile(root,FILE{1}),'TRIAL','GRAND','U','N');

headIdx = 5;
wingIdx = 8;

xIdx = 1;

filename = 'SOS_HeadFree_pat2head_BODE';

FIG = figure (1); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 4 2];
FIG.Name = filename;
movegui(FIG,'center')
hold on

ax1 = subplot(1,1,1) ; hold on
    ax1.FontSize = 8;
    ax1.XLim = [0 10];
    ax1.YLim = [-100 100];
    ax1.XLabel.FontSize = 8;
    ax1.XLabel.Color = 'k';
    ax1.XLabel.String = 'Frequency (Hz)';
 	ax1.YLabel.String = 'Time Difference (ms)';
    ax1.YLabel.FontSize = ax1.XLabel.FontSize;
    
    HFreq           = HeadFree.GRAND{1,headIdx}.Mean{2}{4};
    HPhase          = rad2deg(HeadFree.GRAND{1,headIdx}.CircMean{7}{6}(:,xIdx));
    HPhase_STD      = rad2deg(HeadFree.GRAND{1,headIdx}.CircSTD{7}{6}(:,xIdx));
    HTimeDiff       = 1000*(HPhase/360).*(1./HFreq);
    HTimeDiff_STD   = 1000*(HPhase_STD/360).*(1./HFreq);
    
    [~,h.gain] = PlotPatch(HTimeDiff, HTimeDiff_STD, HFreq, 3, HeadFree.N{1,1}, 'b', [0.4 0.4 0.6], 0.5, 2);
	h.gain.Marker = '.';
    h.gain.MarkerSize = 20;
    
	WFreq           = HeadFree.GRAND{1,wingIdx}.Mean{2}{4};
    WPhase          = rad2deg(HeadFree.GRAND{1,wingIdx}.CircMean{7}{6}(:,xIdx));
    WPhase_STD      = rad2deg(HeadFree.GRAND{1,wingIdx}.CircSTD{7}{6}(:,xIdx));
    WTimeDiff       = 1000*(WPhase/360).*(1./HFreq);
    WTimeDiff_STD   = 1000*(WPhase_STD/360).*(1./HFreq);
    
    [~,h.phase] = PlotPatch(WTimeDiff, WTimeDiff_STD, WFreq, 3, HeadFree.N{1,1}, 'r', [0.4 0.4 0.6], 0.5, 2);
	h.phase.Marker = '.';
    h.phase.MarkerSize = 20;
    
end