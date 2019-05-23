function [FIG] = MakeFig_ChirpLog_HeadFree_pat_head_Pos_Time()
%% MakeFig_ChirpLog_HeadFree_pat_head_Pos_Time.m:
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

HeadFree = load(fullfile(root,FILE{1}),'TRIAL','GRAND','U','N'); % load data structure

figNum = 1;
catIdx = 2;
xIdx = 1;

filename = 'ChirpLog_HeadFree_pat_head_Pos_Time'; % name of figure to save

FIG = figure (figNum); clf % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1200 400]*0.75;
FIG.Name = filename;
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
hold on

% Grand Stats
pp = 1;
for jj = 3
    subplot(1,1,pp) ; hold on ; xlim([0 20])
    title([num2str(HeadFree.U{1,3}{1}(jj)) char(176)],'FontSize',15)
    
    yyaxis right
        ax.R = gca;
        ax.R.YColor = [0 1 0];
        ax.R.XLabel.String = 'Time (s)';
        ax.R.FontSize = 12;
        ylabel(['Stimulus (' char(176) ')'],'FontSize',14)
        ylim(20*[-1 1])
        a = get(ax.R,'XTickLabel');
        set(ax.R,'XTickLabel',a,'FontName','Times','Fontsize',12)

        plot(HeadFree.TRIAL{8,jj}{2,1}.Time,HeadFree.TRIAL{1,jj}{1,1}.X(:,xIdx),'Color','g','LineWidth',2)
    
	yyaxis left
        ax.L = gca;
        ax.L.YColor = [0 0 0];
        ax.L.XLabel.String = 'Time (s)';
        ax.L.FontSize = 12;
        ylabel(['Head (' char(176) ')'],'FontSize',14)
        ylim(20*[-1 1])
        a = get(ax.L,'XTickLabel');
        set(ax.L,'XTickLabel',a,'FontName','Times','Fontsize',12)
        h.patch = PlotPatch(HeadFree.GRAND{jj,catIdx}.Mean{2}{6}(:,xIdx),HeadFree.GRAND{jj,catIdx}.STD{2}{6}(:,xIdx),...
            HeadFree.GRAND{jj,catIdx}.Mean{2}{5},3,HeadFree.N{1,1},'b',[0.4 0.4 0.6],0.5,3);
    
        xlabel('Time (s)','FontSize',14)
        
    pp = pp + 1;
end

%saveas(FIG,[root 'FIGURE\' filename '.fig']); % save .fig file
% print(gcf,[root 'FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
end