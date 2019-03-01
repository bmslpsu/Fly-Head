function [] = MakeFig_Chirp_HeadFree_Head_Pos_Freq(root,figNum)
%% Chirp_HeadFree_Head_Pos_Freq:
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
filename = 'Chirp_HeadFree_Head_Pos_Freq'; % name of figure to save
HeadFree = load([root 'DATA.mat'],'PAT','WING','HEAD','n','unq'); % load data structure

F = figure (figNum); % figure handle
clf
set(gcf,'Position',[200 400 1300 600])
for kk = 1:HeadFree.n.Fly
    for jj = 1:HeadFree.n.Amp
        % MAGNITUDE
        subplot(2,HeadFree.n.Amp,jj) ; hold on
        h.Trial = plot(HeadFree.HEAD.Freq{kk}{jj},HeadFree.HEAD.Mag{kk}{jj},'Color',[0.5 0.5 0.5],'LineWidth',1); % trials
        for ii = 1:length(h.Trial)
            h.Trial(ii).Color(4) = 0.4;
        end
        h.Fly = plot(HeadFree.HEAD.FlyMed.Freq{kk}(:,jj),HeadFree.HEAD.FlyMed.Mag{kk}(:,jj),'LineWidth',1); % flys
        h.Fly.Color(4) = 0.6;
        % PHASE
    	subplot(2,HeadFree.n.Amp,jj+HeadFree.n.Amp) ; hold on
        h.Trial = plot(HeadFree.HEAD.Freq{kk}{jj},HeadFree.HEAD.Phase{kk}{jj},'Color',[0.5 0.5 0.5],'LineWidth',1); % trials
        for ii = 1:length(h.Trial)
            h.Trial(ii).Color(4) = 0.4;
        end
        h.Fly = plot(HeadFree.HEAD.FlyMed.Freq{kk}(:,jj),HeadFree.HEAD.FlyMed.Phase{kk}(:,jj),'LineWidth',1); % flys
        h.Fly.Color(4) = 0.6;
    end
end

for jj = 1:HeadFree.n.Amp
    % MAGNITUDE
	subplot(2,HeadFree.n.Amp,jj) ; hold on
    title([num2str(HeadFree.unq.Amp(jj)) '$^{\circ}$'],'Interpreter','latex','FontSize',15)
        h.patch = PlotPatch(HeadFree.HEAD.GrandMed.Mag(:,jj),HeadFree.HEAD.GrandSTD.Mag(:,jj),...
            HeadFree.HEAD.GrandMed.Freq(:,jj),2,HeadFree.n.Fly,'k',[0.4 0.4 0.6],0.5,2); % all flys
     	xlim([0.1 12])
        ylim([0 4])
        if jj==1
            ylabel('Head Mag($^{\circ}$)','Interpreter','latex','FontSize',15)
        end
        if jj~=1
            yticks(0)
            yticklabels('')
        end
        xticks(0)
        xticklabels('')
    % PHASE    
    subplot(2,HeadFree.n.Amp,jj+HeadFree.n.Amp) ; hold on
    title([num2str(HeadFree.unq.Amp(jj)) '$^{\circ}$'],'Interpreter','latex','FontSize',15)
        h.patch = PlotPatch(HeadFree.HEAD.GrandMed.Phase(:,jj),HeadFree.HEAD.GrandSTD.Phase(:,jj),...
            HeadFree.HEAD.GrandMed.Freq(:,jj),2,HeadFree.n.Fly,'k',[0.4 0.4 0.6],0.5,2); % all flys
     	xlim([0.1 12])
        ylim([-5 5])
        if jj==1
            ylabel('Head Phase(rad)','Interpreter','latex','FontSize',15)
        end
        if jj~=1
            yticks(0)
            yticklabels('')
        end
        xticks([0.1 2:2:12])
        xlabel('Frequency (Hz)','Interpreter','latex','FontSize',15)
end

saveas(F,[root 'FIGURE\' filename '.fig']);
% print('ScreenSizeFigure','-dpng','-r0','-bestfit')
end
