function [FIG] = MakeFig_ChirpLog_HeadFree_pat_head_wing_Time()
%% MakeFig_ChirpLog_HeadFree_pat_head_wing_Time.m:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%

root = 'H:\DATA\Rigid_Data\';

[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

HeadFree = load(fullfile(root,FILE{1}),'TRIAL','GRAND','U','N');

%%
figNum = 1;
catIdx = 2;
xIdx = 1;

filename = 'ChirpLog_HeadFree_pat_head_Time';

FIG = figure (figNum); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 8 2*2];
FIG.Name = filename;
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
hold on
amp = 4;

ax(1) = subplot(2,1,1) ; hold on ; xlim([0 20])
%     title([num2str(HeadFree.U{1,3}{1}(amp)) char(176)],'FontSize',8)
    ax(1).YColor = [0 0 0];
    ax(1).FontSize = 8;
    ax(1).YLabel.String = ['Head (' char(176) ')'];
    ax(1).YLabel.FontSize = 8;
    ax(1).YLim = 20*[-1 1];
    ax(1).YTick = 15*[-1 0 1];
    ax(1).XLabel.String = 'Time (s)';
    ax(1).XLabel.FontSize = 8;
    ax(1).XTick = [ax(1).XLim(1),2:2:ax(1).XLim(2)];

    plot(HeadFree.TRIAL{1,amp}{2,1}.Time,HeadFree.TRIAL{1,amp}{3,1}.X(:,xIdx),'Color','k','LineWidth',1)

    [h.patch] = PlotPatch(HeadFree.GRAND{amp,catIdx}.Mean{1}{6}(:,xIdx),HeadFree.GRAND{amp,catIdx}.STD{1}{6}(:,xIdx),...
        HeadFree.GRAND{amp,catIdx}.Mean{1}{5},1,HeadFree.N{1,1},'b',[0.1 0.1 0.6],0.3,1);
    h.patch.EdgeColor = 'none';

ax(2) = axes; hold on
    ax(2).FontSize = 8;
    ax(2).Position = ax(1).Position + [0 -0.0 0 0];
    ax(2).Color = 'none';
    ax(2).YColor = 'none';
    ax(2).XLabel.String = 'Frequency (Hz)';
    ax(2).XLabel.FontSize = 8;
    ax(2).XAxisLocation = 'top';

    ax(2).XLim = [0.1 12];
    ax(2).XTick = [0.1 0.5 1 2 4 6 8 10 12];
    ax(2).XScale = 'log';

ax(3) = subplot(2,1,2) ; hold on
	ax(3).Color = 'w';
    ax(3).YColor = [0 0 0];
    ax(3).FontSize = 8;
    ax(3).YLabel.String = 'Stimulus (°)';
    ax(3).YAxisLocation = 'right';
    ax(3).YLabel.FontSize = 8;
    ax(3).YLim = 20*[-1 1];
    ax(3).YTick = 20*[-1 0 1];
    ax(3).XLabel.String = 'Time (s)';
    ax(3).XLabel.FontSize = 8;
    ax(3).XTick = [ax(3).XLim(1),2:2:ax(3).XLim(2)];
    
	plot(HeadFree.TRIAL{1,amp}{2,1}.Time,HeadFree.TRIAL{1,amp}{3,1}.X(:,xIdx),'Color','k','LineWidth',1)
    
catIdx = 3;
ax(4) = axes ; hold on ; xlim([0 20])
    title([num2str(HeadFree.U{1,3}{1}(amp)) char(176)],'FontSize',8)
    ax(4).Color = 'none';
    ax(4).YColor = [0 0 0];
    ax(4).FontSize = 8;
    ax(4).YLim = 2*[-1 1];
    ax(4).YTick = 2*[-1 0 1];
	ax(4).YLabel.String = '\DeltaWBA (V)';
    ax(4).YAxisLocation = 'left';
    ax(4).YLabel.FontSize = 8;
    ax(4).XLabel.String = 'Time (s)';
    ax(4).XLabel.FontSize = 8;
    ax(4).XTick = [ax(4).XLim(1),2:2:ax(4).XLim(2)];
    ax(4).Position = ax(3).Position;

    [h.patch] = PlotPatch(HeadFree.GRAND{amp,catIdx}.Mean{1}{6}(:,xIdx),HeadFree.GRAND{amp,catIdx}.STD{1}{6}(:,xIdx),...
        HeadFree.GRAND{amp,catIdx}.Mean{1}{5},1,HeadFree.N{1,1},'r',[0.1 0.1 0.6],0.3,1);
    h.patch.EdgeColor = 'none';
    


end