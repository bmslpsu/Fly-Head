function [] = MakeFig_Chirp_HeadFixed_Wing_Pos_Time(root,figNum)
%% MakeFig_Chirp_HeadFixed_Wing_Pos_Time:
%   INPUTS:
%       root: root directory containing data structure
%       figNum: figure #
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE INPUT %
% root = 'H:\Experiment_HeadExcitation\Chirp\HeadFixed\DATA\';
% figNum = 1;
%---------------------------------------------------------------------------------------------------------------------------------
filename = 'Chirp_HeadFixed_Wing_Pos_Time'; % name of figure to save
HeadFixed = load([root 'Chirp_HeadFixed_DATA.mat'],'PAT','WING','n','unq'); % load data structure

FIG = figure (figNum); % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1100 800];
for kk = 1:HeadFixed.n.Fly
    for jj = 1:HeadFixed.n.Amp
        subplot(HeadFixed.n.Amp,1,jj) ; hold on
        h.Trial = plot(HeadFixed.WING.Time{kk}{jj},HeadFixed.WING.Pos{kk}{jj},'Color',[0.5 0.5 0.5],'LineWidth',1); % trials
        for ii = 1:length(h.Trial)
            h.Trial(ii).Color(4) = 0.25;
        end
        h.Fly = plot(HeadFixed.WING.FlyMed.Time{kk}(:,jj),HeadFixed.WING.FlyMed.Pos{kk}(:,jj),'LineWidth',2); % flys
        h.Fly.Color(4) = 0.5;
    end
end

for jj = 1:HeadFixed.n.Amp
	subplot(HeadFixed.n.Amp,1,jj) ; hold on
    title([num2str(HeadFixed.unq.Amp(jj)) '$^{\circ}$'],'Interpreter','latex','FontSize',15)
        h.patch = PlotPatch(HeadFixed.WING.GrandMed.Pos(:,jj),HeadFixed.WING.GrandSTD.Pos(:,jj),...
            HeadFixed.WING.GrandMed.Time(:,jj),2,HeadFixed.n.Fly,'k',[0.4 0.4 0.6],0.5,3); % all flys
    
     	xlim([0 20])
        ylim(6*[-1 1])
        if jj==HeadFixed.n.Amp
            xlabel('Time (s)','Interpreter','latex','FontSize',15)
        end
        if jj~=HeadFixed.n.Amp
            xticks(0)
            xticklabels('')
        end
        ylabel('$\Delta$WBA(V)','Interpreter','latex','FontSize',15)
end

saveas(FIG,[root 'FIGURE\' filename '.fig']); % save .fig file
print(gcf,[root 'FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
disp('Saved to')
disp(root)
end
