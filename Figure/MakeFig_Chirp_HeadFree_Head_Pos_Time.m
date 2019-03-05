function [] = MakeFig_Chirp_HeadFree_Head_Pos_Time(root,figNum)
%% MakeFig_Chirp_HeadFree_Head_Pos_Time:
%   INPUTS:
%       root: root directory containing data structure
%       figNum: figure #
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
% % EXAMPLE INPUT %
% root = 'H:\Experiment_HeadExcitation\Chirp\HeadFree\DATA\';
% figNum = 1;
%---------------------------------------------------------------------------------------------------------------------------------
filename = 'Chirp_HeadFree_Head_Pos_Time'; % name of figure to save
HeadFree = load([root 'Chirp_HeadFree_DATA.mat'],'PAT','WING','HEAD','n','unq'); % load data structure

FIG = figure (figNum); % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1100 800];
for kk = 1:HeadFree.n.Fly
    for jj = 1:HeadFree.n.Amp
        subplot(HeadFree.n.Amp,1,jj) ; hold on
        h.Trial = plot(HeadFree.HEAD.Time{kk}{jj},HeadFree.HEAD.Pos{kk}{jj},'Color',[0.5 0.5 0.5],'LineWidth',1); % trials
        for ii = 1:length(h.Trial)
            h.Trial(ii).Color(4) = 0.25;
        end
        h.Fly = plot(HeadFree.HEAD.FlyMed.Time{kk}(:,jj),HeadFree.HEAD.FlyMed.Pos{kk}(:,jj),'LineWidth',2); % flys
        h.Fly.Color(4) = 0.5;
    end
end

for jj = 1:HeadFree.n.Amp
	subplot(HeadFree.n.Amp,1,jj) ; hold on
    title([num2str(HeadFree.unq.Amp(jj)) '$^{\circ}$'],'Interpreter','latex','FontSize',15)
        h.patch = PlotPatch(HeadFree.HEAD.GrandMed.Pos(:,jj),HeadFree.HEAD.GrandSTD.Pos(:,jj),...
            HeadFree.HEAD.GrandMed.Time(:,jj),2,HeadFree.n.Fly,'k',[0.4 0.4 0.6],0.5,3); % all flys
    
     	xlim([0 20])
        ylim([-20 20])
        if jj==HeadFree.n.Amp
            xlabel('Time (s)','Interpreter','latex','FontSize',15)
        end
        if jj~=HeadFree.n.Amp
            xticks(0)
            xticklabels('')
        end
        ylabel('Head($^{\circ}$)','Interpreter','latex','FontSize',15)
end

saveas(FIG,[root 'FIGURE\' filename '.fig']); % save .fig file
print(gcf,[root 'FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
disp('Saved to')
disp(root)
end
