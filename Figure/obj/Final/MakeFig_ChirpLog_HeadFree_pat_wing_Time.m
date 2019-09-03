function [FIG] = MakeFig_ChirpLog_HeadFree_pat_wing_Time()
%% MakeFig_ChirpLog_HeadFree_pat_wing_Time:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

HeadFree = load(fullfile(root,FILE{1}),'TRIAL','GRAND','U','N');

clearvars -except HeadFree Amp nAmp

figNum = 1;
catIdx = 3;
xIdx = 1;

filename = 'ChirpLog_HeadFree_pat_wing_Time';

FIG = figure (figNum); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 8 2];
FIG.Name = filename;
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
hold on

pp = 1;
for jj = 3
    ax.R = subplot(1,1,pp) ; hold on ; title([num2str(HeadFree.U{1,3}{1}(jj)) char(176)],'FontSize',15)
    ax.R.YColor = [0 1 0];
    ax.R.XLabel.String = 'Time (s)';
    ax.R.FontSize = 8;
    ax.R.YLabel.String = ['(' char(176) ')'];
    ax.R.YLabel.FontSize = 8;
    ax.R.YLim = 15*[-1 1];
    ax.R.YTick = 15*[-1 0 1];
    ax.R.YAxisLocation = 'right';

	plot(HeadFree.TRIAL{8,jj}{2,1}.Time,HeadFree.TRIAL{10,jj}{1,1}.X(:,xIdx),'Color','g','LineWidth',2)
    
    ax.L = axes; hold on
    ax.L.YColor = [0 0 0];
    ax.L.FontSize = 8;
    ax.L.YLabel.String = '\Delta WBA (V)';
    ax.L.YLabel.FontSize = ax.R.YLabel.FontSize;
    ax.L.YLim = 2*[-1 1];
    ax.L.XLabel.String = 'Time (s)';
    ax.L.XLabel.FontSize = 8;
    ax.L.XColor = 'none';
    ax.L.Color = 'none';
    
    ax.L.Position = ax.R.Position;

    PlotPatch(HeadFree.GRAND{jj,catIdx}.Mean{2}{6}(:,xIdx),HeadFree.GRAND{jj,catIdx}.STD{2}{6}(:,xIdx),...
        HeadFree.GRAND{jj,catIdx}.Mean{2}{5},2,HeadFree.N{1,1},'r',[0.1 0.1 0.6],0.3,2);
     
    pp = pp + 1;
end
end