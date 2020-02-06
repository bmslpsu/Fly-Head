function [FIG] = MakeFig_ChirpLog_HeadFree_pat_head_Time()
%% MakeFig_ChirpLog_HeadFree_pat_head_Time.m:
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
FIG.Position = [1 1 8 2];
FIG.Name = filename;
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
hold on

pp = 1;
for jj = 4
    ax.L = subplot(1,1,pp) ; hold on ; xlim([0 20])
    title([num2str(HeadFree.U{1,3}{1}(jj)) char(176)],'FontSize',8)
    ax.L.YColor = [0 0 0];
    ax.L.FontSize = 8;
    ax.L.YLabel.String = ['Head (' char(176) ')'];
    ax.L.YLabel.FontSize = 8;
    ax.L.YLim = 20*[-1 1];
    ax.L.YTick = 15*[-1 0 1];
    ax.L.XLabel.String = 'Time (s)';
    ax.L.XLabel.FontSize = 8;
    
    plot(HeadFree.TRIAL{1,jj}{2,1}.Time,HeadFree.TRIAL{1,jj}{1,1}.X(:,xIdx),'Color','k','LineWidth',1)

    [h.patch] = PlotPatch(HeadFree.GRAND{jj,catIdx}.Mean{1}{6}(:,xIdx),HeadFree.GRAND{jj,catIdx}.STD{1}{6}(:,xIdx),...
        HeadFree.GRAND{jj,catIdx}.Mean{1}{5},1,HeadFree.N{1,1},'b',[0.1 0.1 0.6],0.3,1);
    h.patch.EdgeColor = 'none';
    
%   	[h.patch] = PlotPatch(HeadFree.GRAND{jj,catIdx}.Median{1}{6}(:,xIdx),HeadFree.GRAND{jj,catIdx}.STD{1}{6}(:,xIdx),...
%         HeadFree.GRAND{jj,catIdx}.Median{1}{5},0,HeadFree.N{1,1},'r',[0.1 0.1 0.6],0.3,1);
%     
%   	[h.patch] = PlotPatch(HeadFree.GRAND{jj,catIdx}.Median{2}{6}(:,xIdx),HeadFree.GRAND{jj,catIdx}.STD{2}{6}(:,xIdx),...
%         HeadFree.GRAND{jj,catIdx}.Median{2}{5},0,HeadFree.N{1,1},'g',[0.1 0.1 0.6],0.3,1);
%     
%     [h.patch] = PlotPatch(HeadFree.GRAND{jj,catIdx}.Mean{2}{6}(:,xIdx),HeadFree.GRAND{jj,catIdx}.STD{2}{6}(:,xIdx),...
%         HeadFree.GRAND{jj,catIdx}.Mean{2}{5},0,HeadFree.N{1,1},'c',[0.1 0.1 0.6],0.3,1);
    
    pp = pp + 1;
end

end