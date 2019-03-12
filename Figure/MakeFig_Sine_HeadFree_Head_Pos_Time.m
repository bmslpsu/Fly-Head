function [] = MakeFig_Sine_HeadFree_Head_Pos_Time(root,figNum)
%% MakeFig_Sine_HeadFree_Head_Pos_Time:
%   INPUTS:
%       root: root directory containing data structure
%       figNum: figure #
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
% % EXAMPLE INPUT %
root = 'H:\EXPERIMENTS\Experiment_Sinusoid\DATA\';
figNum = 1;
%---------------------------------------------------------------------------------------------------------------------------------
filename = 'Sine_HeadFree_Head_Pos_Time'; % name of figure to save
HeadFree{1} = load([root 'Sine_HeadFree_3.75_DATA_.mat'],'PAT','WING','HEAD','BODE','n','unq'); % load data structure
HeadFree{4} = load([root 'Sine_HeadFree_7.5_DATA_.mat'],'PAT','WING','HEAD','BODE','n','unq'); % load data structure
HeadFree{2} = load([root 'Sine_HeadFree_11.25_DATA_.mat'],'PAT','WING','HEAD','BODE','n','unq'); % load data structure
HeadFree{3} = load([root 'Sine_HeadFree_15_DATA_.mat'],'PAT','WING','HEAD','BODE','n','unq'); % load data structure
HeadFree{5} = load([root 'Sine_HeadFree_18.75_DATA_.mat'],'PAT','WING','HEAD','BODE','n','unq'); % load data structure
% nAmp = length(HeadFree);
amp = 5;

FIG = figure (figNum); % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1100 800];
for kk = 1:HeadFree{amp}.n.Fly
    for jj = 1:HeadFree{amp}.n.Freq
        subplot(HeadFree{amp}.n.Freq,1,jj) ; hold on
        h.Trial = plot(HeadFree{amp}.HEAD.Time{kk}{jj},HeadFree{amp}.HEAD.Pos{kk}{jj},'Color',[0.5 0.5 0.5],'LineWidth',1); % trials
        for ii = 1:length(h.Trial)
            h.Trial(ii).Color(4) = 0.25;
        end
        h.Fly = plot(HeadFree{amp}.HEAD.FlyMed.Time{kk}(:,jj),HeadFree{amp}.HEAD.FlyMed.Pos{kk}(:,jj),'LineWidth',2); % flys
        h.Fly.Color(4) = 0.5;
    end
end

for jj = 1:HeadFree{amp}.n.Freq
	subplot(HeadFree{amp}.n.Freq,1,jj) ; hold on
    title([num2str(HeadFree{amp}.unq.Freq(jj)) '$^{\circ}$'],'Interpreter','latex','FontSize',15)
        h.patch = PlotPatch(HeadFree{amp}.HEAD.GrandMed.Pos(:,jj),HeadFree{amp}.HEAD.GrandSTD.Pos(:,jj),...
            HeadFree{amp}.HEAD.GrandMed.Time(:,jj),2,HeadFree{amp}.n.Fly,'k',[0.4 0.4 0.6],0.5,3); % all flys
    
     	xlim([0 10])
        ylim([-20 20])
        if jj==HeadFree{amp}.n.Freq
            xlabel('Time (s)','Interpreter','latex','FontSize',15)
        end
        if jj~=HeadFree{amp}.n.Freq
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
