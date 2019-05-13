function [FIG] = MakeFig_ChirpLog_HeadFree_error_Pos_Time()
%% MakeFig_ChirpLog_HeadFree_error_Pos_Time:
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
catIdx = 4;
xIdx = 1;

filename = 'ChirpLog_HeadFree_error_Pos_Time'; % name of figure to save

FIG = figure (figNum); clf % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1200 800];
FIG.Name = filename;
movegui(FIG,'center')
hold on

% Trials
pp = 1;
for jj = 1:HeadFree.N{1,3}
    for kk = 1:HeadFree.N{1,1}
        for ii = 1:size(HeadFree.TRIAL{kk,jj},1)
            subplot(HeadFree.N{1,3},1,pp) ; hold on ; xlim([0 20]) ; ylim(25*[-1 1])
                h.Trial = plot(HeadFree.TRIAL{kk,jj}{ii,catIdx}.Time,HeadFree.TRIAL{kk,jj}{ii,catIdx}.X(:,xIdx),...
                    '-','Color',[0.5 0.5 0.5 0.1],'LineWidth',1);
        end
    end
    pp = pp + 1;
end

% Fly Stats
pp = 1;
for jj = 1:HeadFree.N{1,3}
    for kk = 1:HeadFree.N{1,1}
        subplot(HeadFree.N{1,3},1,pp) ; hold on ; xlim([0 20]) ; ylim(25*[-1 1])
            h.Fly = plot(HeadFree.FLY{jj}{kk,catIdx}.Mean{5},HeadFree.FLY{jj}{kk,catIdx}.Mean{6}(:,xIdx),...
                '-','LineWidth',1);
            h.Fly.Color(4) = 0.5;
    end
    pp = pp + 1;
end

% Grand Stats
pp = 1;
for jj = 1:HeadFree.N{1,3}
    ax = subplot(HeadFree.N{1,3},1,pp) ; hold on ; xlim([0 20]) ; ylim(25*[-1 1])
    title([num2str(HeadFree.U{1,3}{1}(jj)) char(176)],'FontSize',15)
        h.patch = PlotPatch(HeadFree.GRAND{jj,catIdx}.Mean{2}{6}(:,xIdx),HeadFree.GRAND{jj,catIdx}.STD{2}{6}(:,xIdx),...
            HeadFree.GRAND{jj,catIdx}.Mean{2}{5},3,HeadFree.N{1,1},'k',[0.4 0.4 0.6],0.5,3);
               
    	ylabel('Error ($\circ$)','Interpreter','latex','FontSize',15)

        if pp==HeadFree.N{1,3}
            xlabel('Time (s)','Interpreter','latex','FontSize',15)
       	else
            xticks('')
        end
        
    pp = pp + 1;
end

saveas(FIG,[root 'FIGURE\' filename '.fig']); % save .fig file
% print(gcf,[root 'FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
end