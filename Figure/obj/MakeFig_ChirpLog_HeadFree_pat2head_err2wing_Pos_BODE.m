function [FIG] = MakeFig_ChirpLog_HeadFree_pat2head_err2wing_Pos_BODE()
%% MakeFig_ChirpLog_HeadFree_pat2head_Pos_BODE:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'F:\DATA\Rigid_Data\';

% Select files
[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

HeadFree = load(fullfile(root,FILE{1}),'TRIAL','FLY','GRAND','U','N'); % load data structure
%%
figNum = 1;
catIdx = 5;
xIdx = 1;

filename = 'ChirpLog_HeadFree_pat2head_err2wing_Pos_BODE'; % name of figure to save

FIG = figure (figNum); clf % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1300 600];
FIG.Name = filename;
movegui(FIG,'center')
hold on

% Grand Stats
pp = 1;
for jj = 1:HeadFree.N{1,3}
    ax = subplot(2,HeadFree.N{1,3},pp) ; hold on ; xlim([0.1 7]) ; ylim([0 1.2])
%     title([num2str(HeadFree.U{1,3}{1}(jj)) char(176)],'FontSize',15)
        h.patch = PlotPatch(HeadFree.GRAND{jj,catIdx}.Mean{2}{2}(:,xIdx),HeadFree.GRAND{jj,catIdx}.STD{2}{2}(:,xIdx),...
            HeadFree.GRAND{jj,catIdx}.Mean{2}{1}(:,xIdx),3,HeadFree.N{1,1},'k',[0.4 0.4 0.6],0.5,2);
        
        plot([0 12],[1 1],'--','Color',[0.2 1 0.2 1],'LineWidth',1);
        
        if pp==1
            ylabel(['Head Gain (' char(176) ')'],'FontSize',14)
       	else
            yticks('')
        end
        ax.XTick = sort([0.1 ax.XTick]);
        A_cell = cellfun(@(x) num2str(x),(num2cell(round(2*pi*ax.XTick*HeadFree.U{1,3}{1}(jj)))),'UniformOutput',false);
        
       	a = get(ax,'XTickLabel');
        set(ax,'XTickLabel',a,'FontName','Times','fontsize',12)
        a = get(ax,'YTickLabel');
        set(ax,'YTickLabel',a,'FontName','Times','fontsize',12)
        
        ax2 = axes('Position',get(ax,'Position'),'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k','YColor','k'); 
        k=get(ax,'XTick');
        l=get(ax,'Position');
        set(ax2,'YTick',([]));
        xlim([0.1 12])
      	ax2.XTick = ax.XTick;
        set(ax2,'XTickLabel',A_cell);
        set(ax2,'YTickLabel', num2str(get(ax2,'YTick'),'%.2f\n'),'FontName','Times New Roman','FontSize',12);
        c=get(ax2,'XLabel');
        set(c,'String',['Velocity (' char(176) '/s)']);
        ax.XTick = [];
        
   ax = subplot(2,HeadFree.N{1,3},pp+HeadFree.N{1,3}) ; hold on ; xlim([0.1 7]) ; ylim(pi*[-1 1])
        h.patch = PlotPatch(HeadFree.GRAND{jj,catIdx}.CircMean{9}{3}(:,xIdx),HeadFree.GRAND{jj,catIdx}.CircSTD{9}{3}(:,xIdx),...
            HeadFree.GRAND{jj,catIdx}.Mean{2}{1}(:,xIdx),3,HeadFree.N{1,1},'k',[0.4 0.4 0.6],0.5,2);
        
        plot([0 12],[0 0],'--','Color',[0.2 1 0.2 1],'LineWidth',1);
        
        if pp==1
            ylabel('Head Phase (rad)','FontSize',14)
        else
            yticks('')
        end
        
        ax.XTick = sort([0.1 ax.XTick]);
        xlabel('Frequency (Hz)','FontSize',14)
        
      	a = get(ax,'XTickLabel');
        set(ax,'XTickLabel',a,'FontName','Times','fontsize',12)
        a = get(ax,'YTickLabel');
        set(ax,'YTickLabel',a,'FontName','Times','fontsize',12)

    pp = pp + 1;
end

end