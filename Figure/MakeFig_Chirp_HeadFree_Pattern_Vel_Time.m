function [] = MakeFig_Chirp_HeadFree_Pattern_Vel_Time(root,figNum)
%% MakeFig_Chirp_HeadFree_Pattern_Pos_Time:
%   INPUTS:
%       root: root directory containing data structure
%       figNum: figure #
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE INPUT %
% root = 'E:\Experiment_HeadExcitation\Chirp\HeadFree\DATA\';
% figNum = 1;
%---------------------------------------------------------------------------------------------------------------------------------
filename = 'Chirp_HeadFree_Pattern_Vel_Time'; % name of figure to save
HeadFree = load([root 'DATA.mat'],'PAT','WING','HEAD','n','unq'); % load data structure

F = figure (figNum); % figure handle
clf
set(gcf,'Position',[100 100 1100 800])
for kk = 1:HeadFree.n.Fly
    for jj = 1:HeadFree.n.Amp
        subplot(HeadFree.n.Amp,1,jj) ; hold on
        h.Trial = plot(HeadFree.PAT.Time{kk}{jj},HeadFree.PAT.Vel{kk}{jj},'Color',[0.5 0.5 0.5],'LineWidth',1); % trials
        for ii = 1:length(h.Trial)
            h.Trial(ii).Color(4) = 0.4;
        end
        h.Fly = plot(HeadFree.PAT.FlyMed.Time{kk}(:,jj),HeadFree.PAT.FlyMed.Vel{kk}(:,jj),'LineWidth',1); % flys
        h.Fly.Color(4) = 0.6;
    end
end

for jj = 1:HeadFree.n.Amp
	subplot(HeadFree.n.Amp,1,jj) ; hold on
    title([num2str(HeadFree.unq.Amp(jj)) '$^{\circ}$'],'Interpreter','latex','FontSize',15)
        h.patch = PlotPatch(HeadFree.PAT.GrandMed.Vel(:,jj),HeadFree.PAT.GrandSTD.Vel(:,jj),...
            HeadFree.PAT.GrandMed.Time(:,jj),2,HeadFree.n.Fly,'k',[0.4 0.4 0.6],0.5,2); % all flys
    
     	xlim([0 20])
        ylim(4e4*[-1 1])
        if jj==HeadFree.n.Amp
            xlabel('Time (s)','Interpreter','latex','FontSize',15)
        end
        if jj~=HeadFree.n.Amp
            xticks(0)
            xticklabels('')
        end
        ylabel('Pattern($^{\circ}/s$)','Interpreter','latex','FontSize',15)
end

saveas(F,[root 'FIGURE\' filename '.fig']);
% print('ScreenSizeFigure','-dpng','-r0','-bestfit')
end
