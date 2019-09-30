function [FIG] = MakeFig_SOS_EMD_pat2head_PHASE()
%% MakeFig_SOS_EMD_pat2head_PHASE:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

% Select files
[SOS,SOS_path] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select SOS', root, 'MultiSelect','off');
SOS = cellstr(SOS)';

[EMD,EMD_Path] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select EMD', root, 'MultiSelect','off');
EMD = cellstr(EMD)';

HeadFree    = load(fullfile(SOS_path,SOS{1}),'TRIAL','GRAND','U','N');
EMD         = load(fullfile(EMD_Path,EMD{1}),'Phase','freq');

%%
figNum = 1;
catIdx = 5;
xIdx = 1;

filename = 'SOS_HeadFree_pat2head_BODE';

FIG = figure (figNum); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 4 2];
FIG.Name = filename;
movegui(FIG,'center')
hold on

ax(1) = subplot(1,1,1) ; hold on
    ax(1).FontSize = 8;
    ax(1).XLim = [0 10];
    ax(1).YLim = [-300 60];
   	ax(1).XLabel.String = 'Frequency (Hz)';
    ax(1).XLabel.FontSize = 8;
    ax(1).XLabel.Color = 'k';
 	ax(1).YLabel.String = ['Phase (' char(176) ')'];
    ax(1).YLabel.FontSize = 8;
    ax(1).YTick = -300:60:60;

    [~,h.phase] = PlotPatch(rad2deg(HeadFree.GRAND{1,catIdx}.CircMean{7}{6}(:,xIdx)), rad2deg(HeadFree.GRAND{1,catIdx}.CircSTD{7}{6}(:,xIdx)), ...
        HeadFree.GRAND{1,catIdx}.Mean{2}{4}, 1, HeadFree.N{1,1}, 'k', [0.4 0.4 0.6], 0.5, 2);
    h.phase.Marker = '.';
    h.phase.MarkerSize = 20;

    h.emd = plot(EMD.freq, rad2deg(EMD.Phase(:,4)), 'r', 'LineWidth', 2);
    
    plot(ax(1).XLim,[0 0],'--k','LineWidth',1);
    
    legend([h.phase,h.emd],'SOS Head Phase','EMD Phase')

end