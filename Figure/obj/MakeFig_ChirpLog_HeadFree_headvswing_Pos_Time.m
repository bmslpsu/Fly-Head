function [FIG] = MakeFig_ChirpLog_HeadFree_headvswing_Pos_Time()
%% MakeFig_ChirpLog_HeadFree_headvswing_Pos_Time:
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
catIdx = [2 3];
xIdx = 1;

filename = 'ChirpLog_HeadFree_headvswing_Pos_Time'; % name of figure to save

FIG = figure (figNum); clf % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1400 500];
FIG.Name = filename;
movegui(FIG,'center')
hold on

% Trials
pp = 1;
for jj = 1:HeadFree.N{1,3}
    for kk = 1:HeadFree.N{1,1}
        for ii = 1:size(HeadFree.TRIAL{kk,jj},1)
            subplot(1,HeadFree.N{1,3},pp) ; hold on ; xlim(20*[-1 1]) ; ylim(5*[-1 1])
                h.Trial = plot(HeadFree.TRIAL{kk,jj}{ii,catIdx(1)}.X(:,xIdx),HeadFree.TRIAL{kk,jj}{ii,catIdx(2)}.X(:,xIdx),...
                    '.','Color',[0.5 0.5 0.5 0.1],'LineWidth',1);
        end
    end
    pp = pp + 1;
end

% Fly Stats
pp = 1;
for jj = 1:HeadFree.N{1,3}
    for kk = 1:HeadFree.N{1,1}
        subplot(1,HeadFree.N{1,3},pp) ; hold on ; xlim(20*[-1 1]) ; ylim(5*[-1 1])
            h.Fly = plot(HeadFree.FLY{jj}{kk,catIdx(1)}.Mean{6}(:,xIdx),HeadFree.FLY{jj}{kk,catIdx(2)}.Mean{6}(:,xIdx),...
                '.','LineWidth',1);
            h.Fly.Color(4) = 0.5;
    end
    pp = pp + 1;
end

% Grand Stats
pp = 1;
for jj = 1:HeadFree.N{1,3}
    ax = subplot(1,HeadFree.N{1,3},pp) ; hold on ; xlim(20*[-1 1]) ; ylim(8*[-1 1])
    title([num2str(HeadFree.U{1,3}{1}(jj)) char(176)],'FontSize',15)
        plot(HeadFree.GRAND{jj,catIdx(1)}.Mean{2}{6}(:,xIdx),HeadFree.GRAND{jj,catIdx(2)}.Mean{2}{6}(:,xIdx),'-k','LineWidth',2);
        
        if pp==1
            ylabel('Wing (V)','Interpreter','latex','FontSize',15)
        else
            yticks('')
        end
        
     	xlabel('Head ($\circ$)','Interpreter','latex','FontSize',15)

        
    pp = pp + 1;
end

saveas(FIG,[root 'FIGURE\' filename '.fig']); % save .fig file
% print(gcf,[root 'FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
end