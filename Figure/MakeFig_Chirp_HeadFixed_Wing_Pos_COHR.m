function [] = MakeFig_Chirp_HeadFixed_Wing_Pos_COHR(root,figNum)
%% MakeFig_Chirp_HeadFree_Wing_Pos_COHR:
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
filename = 'Chirp_HeadFixed_Wing_Pos_COHR'; % name of figure to save
HeadFixed = load([root 'Chirp_HeadFixed_DATA.mat'],'PAT','WING','n','unq'); % load data structure

FIG = figure (figNum); % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1100 800];
set(gcf,'Position',[200 400 1300 400])
for kk = 1:HeadFixed.n.Fly
    for jj = 1:HeadFixed.n.Amp
        subplot(1,HeadFixed.n.Amp,jj) ; hold on
        h.Trial = plot(HeadFixed.WING.COHR.Freq{kk}{jj},HeadFixed.WING.COHR.Mag{kk}{jj},'Color',[0.5 0.5 0.5],'LineWidth',1); % trials
        for ii = 1:length(h.Trial)
            h.Trial(ii).Color(4) = 0.2;
        end
        h.Fly = plot(HeadFixed.WING.FlyMean.COHR.Freq{kk}(:,jj),HeadFixed.WING.FlyMean.COHR.Mag{kk}(:,jj),'LineWidth',1); % flys
        h.Fly.Color(4) = 0.6;
    end
end

for jj = 1:HeadFixed.n.Amp
	subplot(1,HeadFixed.n.Amp,jj) ; hold on
    title([num2str(HeadFixed.unq.Amp(jj)) '$^{\circ}$'],'Interpreter','latex','FontSize',15)
        h.patch = PlotPatch(HeadFixed.WING.GrandMean.COHR.Mag(:,jj),HeadFixed.WING.GrandSTD.COHR.Mag(:,jj),...
            HeadFixed.WING.GrandMean.COHR.Freq(:,jj),2,HeadFixed.n.Fly,'k',[0.4 0.4 0.6],0.5,2); % all flys
        xlim([0.1 12])
        ylim([0 1])
        if jj==1
            ylabel('Wing Coherence','Interpreter','latex','FontSize',15)
        end
        if jj~=1
            yticks(0)
            yticklabels('')
        end
        xticks([0.1 2:2:12])
        xlabel('Frequency (Hz)','Interpreter','latex','FontSize',15)
end

saveas(FIG,[root 'FIGURE\' filename '.fig']); % save .fig file
print(gcf,[root 'FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
disp('Saved to')
disp(root)
end
