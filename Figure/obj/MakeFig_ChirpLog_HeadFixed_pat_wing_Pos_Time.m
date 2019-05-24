function [FIG] = MakeFig_ChirpLog_HeadFixed_pat_wing_Pos_Time()
%% MakeFig_ChirpLog_HeadFixed_pat_wing_Pos_Time.m:
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

HeadFixed = load(fullfile(root,FILE{1}),'TRIAL','GRAND','U','N'); % load data structure
%%
figNum = 1;
catIdx = 2;
xIdx = 1;

filename = 'ChirpLog_HeadFixed_pat_wing_Pos_Time'; % name of figure to save

FIG = figure (figNum); clf % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 900 600/2];
FIG.Name = filename;
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
hold on

% Grand Stats
for jj = 3
	hold on ; xlim([0 20])
    title([num2str(HeadFixed.U{1,3}{1}(jj)) char(176)],'FontSize',15)
    
    yyaxis right
        ax.R = gca;
        ax.R.YColor = [0 1 0];
        ax.R.XLabel.String = 'Time (s)';
        ax.R.FontSize = 12;
        ylabel(['Stimulus (' char(176) ')'],'FontSize',14)
        ylim(20*[-1 1])
        a = get(ax.R,'XTickLabel');
        set(ax.R,'XTickLabel',a,'FontName','Times','fontsize',12)
        a = get(ax.R,'YTickLabel');
        set(ax.R,'YTickLabel',a,'FontName','Times','fontsize',12)
        plot(HeadFixed.TRIAL{8,jj}{2,1}.Time,HeadFixed.TRIAL{1,jj}{1,1}.X(:,xIdx),'Color','g','LineWidth',2)
        
%         h.patch = PlotPatch(HeadFree.GRAND{jj,2}.Mean{2}{6}(:,xIdx),HeadFree.GRAND{jj,2}.STD{2}{6}(:,xIdx),...
%             HeadFree.GRAND{jj,2}.Mean{2}{5},2,HeadFree.N{1,1},'b',[0.4 0.4 0.6],0.5,3);
        
	yyaxis left
        ax.L = gca;
        ax.L.YColor = [0 0 0];
        ax.L.XLabel.String = 'Time (s)';
        ax.L.FontSize = 12;
        ylabel('\Delta WBA (V)','FontSize',14)
        ylim(3*[-1 1])
        a = get(ax.L,'XTickLabel');
        set(ax.L,'XTickLabel',a,'FontName','Times','fontsize',12)
       	a = get(ax.L,'YTickLabel');
        set(ax.L,'YTickLabel',a,'FontName','Times','fontsize',12)
        h.patch = PlotPatch(HeadFixed.GRAND{jj,catIdx}.Mean{2}{6}(:,xIdx),HeadFixed.GRAND{jj,catIdx}.STD{2}{6}(:,xIdx),...
            HeadFixed.GRAND{jj,catIdx}.Mean{2}{5},3,HeadFixed.N{1,1},'r',[0.4 0.4 0.6],0.5,3);
        
        xlabel('Time (s)','FontSize',14)
end

end