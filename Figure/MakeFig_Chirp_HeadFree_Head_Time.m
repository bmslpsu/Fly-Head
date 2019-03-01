function [] = MakeFig_Chirp_HeadFree_Head_Time(root,figNum)
%% Fig_Chirp_HeadFree_Time:
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
filename = 'MakeFig_Chirp_HeadFree_Head_Time'; % name of file to save
load([root 'DATA.mat'],'PAT','WING','HEAD','n','unq') % load data structure

F = figure (figNum); % figure handle
clf
set(gcf,'Position',[100 100 1100 800])
for kk = 1:n.Fly
    for jj = 1:n.Amp
        subplot(n.Amp,1,jj) ; hold on
        h.Trial = plot(HEAD.Time{kk}{jj},HEAD.Pos{kk}{jj},'Color',[0.5 0.5 0.5],'LineWidth',1); % all trials
        for ii = 1:length(h.Trial)
            h.Trial(ii).Color(4) = 0.25;
        end
        h.Fly = plot(HEAD.FlyMed.Time{kk}(:,jj),HEAD.FlyMed.Pos{kk}(:,jj),'LineWidth',2); % all trials
        h.Fly.Color(4) = 0.5;
    end
end

for jj = 1:n.Amp
	subplot(n.Amp,1,jj) ; hold on ; title([num2str(unq.Amp(jj)) '$^{\circ}$'],'Interpreter','latex','FontSize',15)
        h.patch = PlotPatch(HEAD.GrandMed.Pos(:,jj),HEAD.GrandSTD.Pos(:,jj),HEAD.GrandMed.Time(:,jj),2,n.Fly,'k',...
            [0.4 0.4 0.6],0.7,3);
    
     	xlim([0 20])
        ylim([-20 20])
        if jj==n.Amp
            xlabel('Time (s)','Interpreter','latex','FontSize',15)
        end
        if jj~=n.Amp
            xticks(0)
            xticklabels('')
        end
        ylabel('Head($^{\circ}$)','Interpreter','latex','FontSize',15)
end

saveas(F,[root 'FIGURE\' filename '.fig']);
% print('ScreenSizeFigure','-dpng','-r0','-bestfit')
end
