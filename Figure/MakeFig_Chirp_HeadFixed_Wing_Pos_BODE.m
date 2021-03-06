function [] = MakeFig_Chirp_HeadFixed_Wing_Pos_BODE(root,figNum)
%% MakeFig_Chirp_HeadFixed_Wing_Pos_BODE:
%   INPUTS:
%       root: root directory containing data structure
%       figNum: figure #
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE INPUT %
root = 'H:\Experiment_HeadExcitation\Chirp\HeadFixed\DATA\';
figNum = 1;
%---------------------------------------------------------------------------------------------------------------------------------
filename = 'Chirp_HeadFixed_Wing_Pos_BODE'; % name of figure to save
HeadFixed = load([root 'Chirp_HeadFixed_DATA.mat'],'PAT','WING','n','unq'); % load data structure

FIG = figure (figNum); % figure handle
FIG.Color = 'w';
set(gcf,'Position',[200 400 1300 600])
for kk = 1:HeadFixed.n.Fly
    for jj = 1:HeadFixed.n.Amp
        % MAGNITUDE
        subplot(2,HeadFixed.n.Amp,jj) ; hold on
        h.Trial = plot(HeadFixed.WING.Freq{kk}{jj},HeadFixed.WING.GAIN{kk}{jj},'Color',[0.5 0.5 0.5],'LineWidth',1); % trials
        for ii = 1:length(h.Trial)
            h.Trial(ii).Color(4) = 0.2;
        end
        h.Fly = plot(HeadFixed.WING.FlyMean.Freq{kk}(:,jj),HeadFixed.WING.FlyMean.GAIN{kk}(:,jj),'LineWidth',1); % flys
        h.Fly.Color(4) = 0.6;
        % PHASE
    	subplot(2,HeadFixed.n.Amp,jj+HeadFixed.n.Amp) ; hold on
        h.Trial = plot(HeadFixed.WING.Freq{kk}{jj},HeadFixed.WING.PHASE{kk}{jj},'Color',[0.5 0.5 0.5],'LineWidth',1); % trials
        for ii = 1:length(h.Trial)
            h.Trial(ii).Color(4) = 0.4;
        end
        h.Fly = plot(HeadFixed.WING.FlyMean.Freq{kk}(:,jj),HeadFixed.WING.FlyMean.PHASE{kk}(:,jj),'LineWidth',1); % flys
        h.Fly.Color(4) = 0.6;
    end
end

for jj = 1:HeadFixed.n.Amp
    % MAGNITUDE
	subplot(2,HeadFixed.n.Amp,jj) ; hold on
    title([num2str(HeadFixed.unq.Amp(jj)) '$^{\circ}$'],'Interpreter','latex','FontSize',15)
        h.patch = PlotPatch(HeadFixed.WING.GrandMean.GAIN(:,jj),HeadFixed.WING.GrandSTD.GAIN(:,jj),...
            HeadFixed.WING.GrandMean.Freq(:,jj),2,HeadFixed.n.Fly,'k',[0.4 0.4 0.6],0.5,2); % all flys
     	xlim([0.1 12])
        ylim([0 0.5])
        if jj==1
            ylabel('Wing Gain($V/{\circ}$)','Interpreter','latex','FontSize',15)
        end
        if jj~=1
            yticks(0)
            yticklabels('')
        end
        xticks(0)
        xticklabels('')
    % PHASE    
    subplot(2,HeadFixed.n.Amp,jj+HeadFixed.n.Amp) ; hold on
        h.patch = PlotPatch(HeadFixed.WING.GrandMean.PHASE(:,jj),HeadFixed.WING.GrandSTD.PHASE(:,jj),...
            HeadFixed.WING.GrandMean.Freq(:,jj),2,HeadFixed.n.Fly,'k',[0.4 0.4 0.6],0.5,2); % all flys
     	xlim([0.1 12])
        ylim([-4 4])
        if jj==1
            ylabel('Wing Phase(rad)','Interpreter','latex','FontSize',15)
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
