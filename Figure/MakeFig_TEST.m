function [] = MakeFig_TEST(root,figNum)
%% MakeFig_Chirp_HeadFree_Pattern_Pos_Time:
%   INPUTS:
%       root: root directory containing data structure
%       figNum: figure #
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE INPUT %
root = 'H:\EXPERIMENTS\Experiment_ChirpLog_HeadFree\DATA\';
figNum = 1;
%---------------------------------------------------------------------------------------------------------------------------------
filename = 'Chirp_HeadFree_Pattern_Pos_Time'; % name of figure to save
HeadFree = load([root 'Chirp_HeadFree_DATA_TEST.mat'],'TRIAL','FLY','GRAND','U','N'); % load data structure
%%
FIG = figure (figNum); % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1100 800];

% Trials
for kk = 1:HeadFree.N{1,1}
    for jj = 1:HeadFree.N{1,3}
        for ii = size(HeadFree.TRIAL{kk}{jj},1)
            subplot(HeadFree.N{1,3},1,jj) ; hold on
            h.Trial = plot(HeadFree.TRIAL{kk}{jj}{ii,2}.Time,HeadFree.TRIAL{kk}{jj}{ii,2}.X(:,1),'Color',[0.5 0.5 0.5],'LineWidth',1); % trials
         	h.Trial.Color(4) = 0.25;
        end
    end
end

% Fly Stats
for kk = 1:HeadFree.N{1,1}
    for jj = 1:HeadFree.N{1,3}
        subplot(HeadFree.N{1,3},1,jj) ; hold on
        h.Fly = plot(HeadFree.FLY{kk}{jj}.Mean{5,2},HeadFree.FLY{kk}{jj}.Mean{6,2}(:,1),'LineWidth',2); % flys
        h.Fly.Color(4) = 0.5;
    end
end

% Grand Stats
for jj = 1:HeadFree.N{1,3}
	subplot(HeadFree.N{1,3},1,jj) ; hold on
    title([num2str(HeadFree.U{1,3}{1}(jj)) '$^{\circ}$'],'Interpreter','latex','FontSize',15)
        h.patch = PlotPatch(HeadFree.GRAND{jj}.Mean{2}{6,2}(:,1),HeadFree.GRAND{jj}.STD{2}{6,2}(:,1),...
            HeadFree.GRAND{jj}.Mean{2}{5,2},2,HeadFree.N{1,1},'k',[0.4 0.4 0.6],0.5,3); % all flys
    
     	xlim([0 20])
        ylim([-20 20])
        if jj==HeadFree.N{1,3}
            xlabel('Time (s)','Interpreter','latex','FontSize',15)
        end
        if jj~=HeadFree.N{1,3}
            xticks(0)
            xticklabels('')
        end
        ylabel('Head($^{\circ}$)','Interpreter','latex','FontSize',15)
end

% saveas(FIG,[root 'FIGURE\' filename '.fig']); % save .fig file
% print(gcf,[root 'FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
% disp('Saved to')
% disp(root)
end