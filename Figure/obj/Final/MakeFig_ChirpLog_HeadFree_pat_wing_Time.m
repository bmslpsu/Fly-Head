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

figNum = 1;
catIdx = 3;
xIdx = 1;

filename = 'ChirpLog_HeadFree_pat_wing_Time';

FIG = figure (figNum); clf
FIG.Color = 'w';
FIG.Position = [100 100 1200 400]*0.75;
FIG.Name = filename;
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
hold on

pp = 1;
for jj = 3
    subplot(1,1,pp) ; hold on ; xlim([0 20])
    title([num2str(HeadFree.U{1,3}{1}(jj)) char(176)],'FontSize',15)
    
    yyaxis right
        ax.R = gca;
        ax.R.YColor = [0 1 0];
        ax.R.XLabel.String = 'Time (s)';
        ax.R.FontSize = 12;
        ax.R.YLabel.String = ['Stimulus (' char(176) ')'];
        ax.R.YLabel.FontSize = 14;
        ax.R.YLim = 20*[-1 1];

        plot(HeadFree.TRIAL{8,jj}{2,1}.Time,HeadFree.TRIAL{10,jj}{1,1}.X(:,xIdx),'Color','g','LineWidth',2)
    
	yyaxis left
        ax.L = gca;
        ax.L.YColor = [0 0 0];
        ax.L.FontSize = ax.R.FontSize;
     	ax.L.YLabel.String = ['Wings (' char(176) ')'];
        ax.L.YLabel.FontSize = ax.R.YLabel.FontSize;
        ax.L.YLim = 2*[-1 1];
     	ax.L.XLabel.String = 'Time (s)';
        ax.L.XLabel.FontSize = ax.R.YLabel.FontSize;

        PlotPatch(HeadFree.GRAND{jj,catIdx}.Mean{2}{6}(:,xIdx),HeadFree.GRAND{jj,catIdx}.STD{2}{6}(:,xIdx),...
            HeadFree.GRAND{jj,catIdx}.Mean{2}{5},3,HeadFree.N{1,1},'r',[0.1 0.1 0.6],0.3,3);
     
    pp = pp + 1;
end
end