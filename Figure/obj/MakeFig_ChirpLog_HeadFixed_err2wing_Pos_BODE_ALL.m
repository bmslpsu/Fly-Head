function [FIG] = MakeFig_ChirpLog_HeadFixed_err2wing_Pos_BODE_ALL()
%% MakeFig_ChirpLog_HeadFixed_err2wing_Pos_BODE_ALL:
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

HeadFree = load(fullfile(root,FILE{1}),'TRIAL','FLY','GRAND','U','N'); % load data structure

figNum = 1;
catIdx = 3;
xIdx = 1;

filename = 'ChirpLog_HeadFixed_err2wing_Pos_BODE_ALL'; % name of figure to save

FIG = figure (figNum); clf % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1100 700];
FIG.Name = filename;
movegui(FIG,'center')
hold on
ax = gca;

CLR         = jet(HeadFree.N{1,3});
legList     = cellstr(num2str(HeadFree.U{1,3}{1}))';
hh          = cell(HeadFree.N{1,3},1);

% Grand Stats
pp = 1;
for jj = 1:HeadFree.N{1,3}
    ax = subplot(2,1,1) ; hold on ; xlim([0.1 12]) ; ylim([0 0.6])
    title('All','FontSize',15)
        [~,h.patch] = PlotPatch(HeadFree.GRAND{jj,catIdx}.Mean{2}{2}(:,xIdx),HeadFree.GRAND{jj,catIdx}.STD{2}{2}(:,xIdx),...
            HeadFree.GRAND{jj,catIdx}.Mean{2}{1}(:,xIdx),3,HeadFree.N{1,1},CLR(jj,:),[0.4 0.4 0.6],0.5,3);
        hh{jj} = h.patch;
        plot([0 12],[1 1],'--','Color',[0.2 1 0.2 1],'LineWidth',1);
        
        if pp==1
            ylabel('Wing Gain ($V/^{\circ}$)','Interpreter','latex','FontSize',15)
       	else
%             yticks('')
        end
        
        xticks('')

   ax = subplot(2,1,2) ; hold on ; xlim([0.1 12]) ; ylim(pi*[-1 1])
        h.patch = PlotPatch(HeadFree.GRAND{jj,catIdx}.CircMean{9}{3}(:,xIdx),HeadFree.GRAND{jj,catIdx}.CircSTD{9}{3}(:,xIdx),...
            HeadFree.GRAND{jj,catIdx}.Mean{2}{1}(:,xIdx),2,HeadFree.N{1,1},CLR(jj,:),[0.4 0.4 0.6],0.5,3);
        
        plot([0 12],[0 0],'--','Color',[0.2 1 0.2 1],'LineWidth',1);
        
        if pp==1
            ylabel('Head Phase (rad)','Interpreter','latex','FontSize',15)
        else
%             yticks('')
        end
        
        ax.XTick = unique(sort([0.1 ax.XTick]));
        xlabel('Frequency (Hz)','Interpreter','latex','FontSize',15)

    pp = pp + 1;
end
subplot(2,1,1) ; hold on
leg = legend(cat(1,hh{:}),legList);
leg.Title.String = ['Amplitude' char(176)];

saveas(FIG,[root 'FIGURE\' filename '.fig']); % save .fig file
% print(gcf,[root 'FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
end