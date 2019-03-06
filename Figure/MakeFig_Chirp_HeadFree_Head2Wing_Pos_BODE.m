function [] = MakeFig_Chirp_HeadFree_Head2Wing_Pos_BODE(root,figNum)
%% MakeFig_Chirp_HeadFree_Head2Wing_Pos_BODE:
%   INPUTS:
%       root: root directory containing data structure
%       figNum: figure #
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE INPUT %
% root = 'E:\EXPERIMENTS\Experiment_ChirpLog_HeadFree\DATA\';
% figNum = 1;
%---------------------------------------------------------------------------------------------------------------------------------
filename = 'Chirp_HeadFree_Head2Wing_Pos_BODE'; % name of figure to save
HeadFree = load([root 'Chirp_HeadFree_DATA.mat'],'PAT','WING','HEAD','BODE','n','unq'); % load data structure

FIG = figure (figNum); % figure handle
FIG.Color = 'w';
set(gcf,'Position',[200 400 1300 600])
% for kk = 1:HeadFree.n.Fly
%     for jj = 1:HeadFree.n.Amp
%         % MAGNITUDE
%         subplot(2,HeadFree.n.Amp,jj) ; hold on
% %         h.Trial = plot(HeadFree.HEAD.Freq{kk}{jj},HeadFree.BODE.head2wing.GAIN{kk}{jj},'Color',[0.5 0.5 0.5],'LineWidth',1); % trials
% %         for ii = 1:length(h.Trial)
% %             h.Trial(ii).Color(4) = 0.2;
% %         end
%         h.Fly = plot(HeadFree.HEAD.FlyMean.Freq{kk}(:,jj),HeadFree.BODE.FlyMean.head2wing.GAIN{kk}(:,jj),'LineWidth',1); % flys
%         h.Fly.Color(4) = 0.6;
%         % PHASE
%     	subplot(2,HeadFree.n.Amp,jj+HeadFree.n.Amp) ; hold on
% %         h.Trial = plot(HeadFree.HEAD.Freq{kk}{jj},HeadFree.BODE.head2wing.PHASE{kk}{jj},'Color',[0.5 0.5 0.5],'LineWidth',1); % trials
% %         for ii = 1:length(h.Trial)
% %             h.Trial(ii).Color(4) = 0.4;
% %         end
%         h.Fly = plot(HeadFree.HEAD.FlyMean.Freq{kk}(:,jj),HeadFree.BODE.FlyMean.head2wing.PHASE{kk}(:,jj),'LineWidth',1); % flys
%         h.Fly.Color(4) = 0.6;
%     end
% end

for jj = 1:HeadFree.n.Amp
    % MAGNITUDE
	subplot(2,HeadFree.n.Amp,jj) ; hold on
    title([num2str(HeadFree.unq.Amp(jj)) '$^{\circ}$'],'Interpreter','latex','FontSize',15)
        h.patch = PlotPatch(HeadFree.BODE.GrandMean.head2wing.GAIN(:,jj),HeadFree.BODE.GrandSTD.head2wing.GAIN(:,jj),...
            HeadFree.HEAD.GrandMean.Freq(:,jj),2,HeadFree.n.Fly,'k',[0.4 0.4 0.6],0.5,2); % all flys
     	xlim([0.1 12])
        ylim([0 1])
        if jj==1
            ylabel('Head Gain(${\circ}/{\circ}$)','Interpreter','latex','FontSize',15)
        end
        if jj~=1
            yticks(0)
            yticklabels('')
        end
        xticks(0)
        xticklabels('')
    % PHASE    
    subplot(2,HeadFree.n.Amp,jj+HeadFree.n.Amp) ; hold on
        h.patch = PlotPatch(HeadFree.BODE.GrandMean.head2wing.PHASE(:,jj),HeadFree.BODE.GrandSTD.head2wing.PHASE(:,jj),...
            HeadFree.HEAD.GrandMean.Freq(:,jj),2,HeadFree.n.Fly,'k',[0.4 0.4 0.6],0.5,2); % all flys
     	xlim([0.1 12])
        ylim([-2 2])
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

saveas(FIG,[root 'FIGURE\' filename '.fig']); % save .fig file
print(gcf,[root 'FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
disp('Saved to')
disp(root)
end
