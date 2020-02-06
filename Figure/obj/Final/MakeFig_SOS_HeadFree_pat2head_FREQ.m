function [FIG] = MakeFig_SOS_HeadFree_pat2head_FREQ()
%% MakeFig_SOS_HeadFree_pat2head_FREQ:
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
catIdx = 2;
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
    ax1.YLim = [0 10];
    ax1.XTickLabel = '';
    ax1.XLabel.FontSize = 8;
 	ax1.YLabel.String = ['Position (' char(176) ')'];
    ax1.YLabel.FontSize = ax1.XLabel.FontSize;
    
    h.ref = plot(HeadFree.TRIAL{1}{1,1}.Fv, HeadFree.TRIAL{1}{1,1}.Mag(:,xIdx), 'k','LineWidth',1);
    
    [~,h.mag] = PlotPatch(HeadFree.GRAND{1,catIdx}.Mean{2}{8}(:,xIdx), HeadFree.GRAND{1,catIdx}.STD{2}{8}(:,xIdx), ...
        HeadFree.GRAND{1,catIdx}.Mean{2}{7}, 1, HeadFree.N{1,1}, 'b', [0.4 0.4 0.6], 0.5, 1);
    
    leg = legend([h.ref,h.mag],'Stimulus','Head');
    leg.Box = 'off';

ax2 = subplot(2,1,2) ; hold on
    ax2.FontSize = ax1.FontSize;
    ax2.XLim = ax1.XLim;
    ax2.YLim = [0 100];
   	ax2.XLabel.String = 'Frequency (Hz)';
    ax2.XLabel.FontSize = 8;
    ax2.XLabel.Color = 'k';
 	ax2.YLabel.String = ['Velocity (' char(176) '/s)'];
    ax2.YLabel.FontSize = ax1.XLabel.FontSize;

    plot(HeadFree.TRIAL{1}{1,1}.Fv, HeadFree.TRIAL{1}{1,1}.Mag(:,xIdx+1), 'k','LineWidth',1);
    
    [~,h.mag] = PlotPatch(HeadFree.GRAND{1,catIdx}.Mean{2}{8}(:,xIdx+1), HeadFree.GRAND{1,catIdx}.STD{2}{8}(:,xIdx+1), ...
        HeadFree.GRAND{1,catIdx}.Mean{2}{7}, 1, HeadFree.N{1,1}, 'b', [0.4 0.4 0.6], 0.5, 1);

end