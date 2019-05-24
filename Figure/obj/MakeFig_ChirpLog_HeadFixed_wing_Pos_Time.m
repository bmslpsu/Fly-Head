function [FIG] = MakeFig_ChirpLog_HeadFixed_wing_Pos_Time()
%% MakeFig_ChirpLog_HeadFixed_wing_Pos_Time:
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

HeadFixed = load(fullfile(root,FILE{1}),'TRIAL','FLY','GRAND','U','N'); % load data structure

figNum = 1;
catIdx = 2;
xIdx = 1;

filename = 'ChirpLog_HeadFixed_wing_Pos_Time'; % name of figure to save

FIG = figure (figNum); clf % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 900 600/4];
FIG.Name = filename;
movegui(FIG,'center')
hold on

% % Trials
% pp = 1;
% for jj = 1:HeadFixed.N{1,3}
%     for kk = 1:HeadFixed.N{1,1}
%         for ii = 1:size(HeadFixed.TRIAL{kk,jj},1)
%             subplot(HeadFixed.N{1,3},1,pp) ; hold on ; xlim([0 20]) ; ylim(3*[-1 1])
%                 h.Trial = plot(HeadFixed.TRIAL{kk,jj}{ii,catIdx}.Time,HeadFixed.TRIAL{kk,jj}{ii,catIdx}.X(:,xIdx),...
%                     '-','Color',[0.5 0.5 0.5 0.1],'LineWidth',1);
%         end
%     end
%     pp = pp + 1;
% end
% 
% % Fly Stats
% pp = 1;
% for jj = 1:HeadFixed.N{1,3}
%     for kk = 1:HeadFixed.N{1,1}
%         subplot(HeadFixed.N{1,3},1,pp) ; hold on ; xlim([0 20]) ; ylim(3*[-1 1])
%             h.Fly = plot(HeadFixed.FLY{jj}{kk,catIdx}.Mean{5},HeadFixed.FLY{jj}{kk,catIdx}.Mean{6}(:,xIdx),...
%                 '-','LineWidth',1);
%             h.Fly.Color(4) = 0.5;
%     end
%     pp = pp + 1;
% end

% Grand Stats
pp = 1;
for jj = 3
    ax = subplot(1,1,pp) ; hold on ; xlim([0 20]) ; ylim(4*[-1 1])
    title([num2str(HeadFixed.U{1,3}{1}(jj)) char(176)],'FontSize',14)
        h.patch = PlotPatch(HeadFixed.GRAND{jj,catIdx}.Mean{2}{6}(:,xIdx),HeadFixed.GRAND{jj,catIdx}.STD{2}{6}(:,xIdx),...
            HeadFixed.GRAND{jj,catIdx}.Mean{2}{5},3,HeadFixed.N{1,1},'c',[0.4 0.4 0.6],0.5,3);
               
    	ylabel('\Delta WBA (V)','FontSize',14)
        ax.FontSize = 12;

        if pp==HeadFixed.N{1,3}
            xlabel('Time (s)','FontSize',14)
       	else
            ax.XTickLabels = '';
        end
        
    pp = pp + 1;
end

end