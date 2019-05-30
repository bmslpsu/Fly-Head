function [FIG] = MakeFig_ChirpLog_HeadFixed_err2wing_Pos_BODE()
%% MakeFig_ChirpLog_HeadFixed_err2wing_Pos_BODE:
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

figNum = 1;
catIdx = 3;
xIdx = 1;

filename = 'ChirpLog_HeadFixed_err2wing_Pos_BODE'; % name of figure to save

FIG = figure (figNum); clf % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1100 700];
FIG.Name = filename;
movegui(FIG,'center')
hold on

% % Trials
% pp = 1;
% for jj = 1:HeadFree.N{1,3}
%     for kk = 1:HeadFree.N{1,1}
%         for ii = 1:size(HeadFree.TRIAL{kk,jj},1)
%             subplot(2,HeadFree.N{1,3},pp) ; hold on ; xlim([0 12.]) ; ylim([0 0.5])
%                 h.Trial = plot(HeadFree.TRIAL{kk,jj}{ii,catIdx}.BodeFv(:,xIdx),HeadFree.TRIAL{kk,jj}{ii,catIdx}.BodeGain(:,xIdx),...
%                     '-','Color',[0.5 0.5 0.5 0.1],'LineWidth',1);
%             subplot(2,HeadFree.N{1,3},pp+HeadFree.N{1,3}) ; hold on ; xlim([0 12.]) ; ylim(pi*[-1 1])
%                 h.Trial = plot(HeadFree.TRIAL{kk,jj}{ii,catIdx}.BodeFv(:,xIdx),HeadFree.TRIAL{kk,jj}{ii,catIdx}.BodePhaseDiff(:,xIdx),...
%                     '-','Color',[0.5 0.5 0.5 0.1],'LineWidth',1);
%         end
%     end
%     pp = pp + 1;
% end


% % Fly Stats
% pp = 1;
% for jj = 1:HeadFree.N{1,3}
%     for kk = 1:HeadFree.N{1,1}
%         subplot(2,HeadFree.N{1,3},pp) ; hold on ; xlim([0 12.]) ; ylim([0 0.5])
%             h.Fly = plot(HeadFree.FLY{jj}{kk,catIdx}.Mean{1}(:,xIdx),HeadFree.FLY{jj}{kk,catIdx}.Median{2}(:,xIdx),...
%                 '-','LineWidth',1);
%             h.Fly.Color(4) = 0.5;
%         subplot(2,HeadFree.N{1,3},pp+HeadFree.N{1,3}) ; hold on ; xlim([0 12.]) ; ylim(pi*[-1 1])
%             h.Fly = plot(HeadFree.FLY{jj}{kk,catIdx}.Mean{1}(:,xIdx),HeadFree.FLY{jj}{kk,catIdx}.CircMean{3}(:,xIdx),...
%                 '-','LineWidth',1);
%             h.Fly.Color(4) = 0.5;
%     end
%     pp = pp + 1;
% end

% Grand Stats
pp = 1;
for jj = 1:HeadFree.N{1,3}
    ax = subplot(2,HeadFree.N{1,3},pp) ; hold on ; xlim([0.1 12]) ; ylim([0 0.6])
    title([num2str(HeadFree.U{1,3}{1}(jj)) char(176)],'FontSize',15)
        h.patch = PlotPatch(medfilt1(HeadFree.GRAND{jj,catIdx}.Mean{2}{2}(:,xIdx),1),medfilt1(HeadFree.GRAND{jj,catIdx}.STD{2}{2}(:,xIdx),1),...
            HeadFree.GRAND{jj,catIdx}.Mean{2}{1}(:,xIdx),3,HeadFree.N{1,1},'k',[0.4 0.4 0.6],0.5,2);
               
        if pp==1
            ylabel('Wing Gain ($V/^{\circ}$)','Interpreter','latex','FontSize',15)
       	else
            yticks('')
        end
        
        xticks('')

    ax = subplot(2,HeadFree.N{1,3},pp+HeadFree.N{1,3}) ; hold on ; xlim([0.1 12]) ; ylim(pi*[-1 1])
        h.patch = PlotPatch(medfilt1(HeadFree.GRAND{jj,catIdx}.CircMean{9}{3}(:,xIdx),3),medfilt1(HeadFree.GRAND{jj,catIdx}.CircSTD{9}{3}(:,xIdx),5),...
            HeadFree.GRAND{jj,catIdx}.Mean{2}{1}(:,xIdx),3,HeadFree.N{1,1},'k',[0.4 0.4 0.6],0.5,2);
        
        plot([0 12],[0 0],'--','Color',[0.2 1 0.2 1],'LineWidth',1);
        
        if pp==1
            ylabel('Wing Phase (rad)','Interpreter','latex','FontSize',15)
        else
            yticks('')
        end
        ax.XTick = sort([0.1 ax.XTick]);
        xlabel('Frequency (Hz)','Interpreter','latex','FontSize',15)

    pp = pp + 1;
end

saveas(FIG,[root 'FIGURE\' filename '.fig']); % save .fig file
% print(gcf,[root 'FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
end