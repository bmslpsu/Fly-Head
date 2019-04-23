function [] = MakeFig_Chirp_Walking_Head_Pos_Time(root,figNum)
%% MakeFig_Chirp_HeadFree_Head_Pos_Time:
%   INPUTS:
%       root: root directory containing data structure
%       figNum: figure #
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE INPUT %
root = 'W:\Research\Walking Chirp mat';
figNum = 1;
%---------------------------------------------------------------------------------------------------------------------------------
filename = 'Chirp_Walking_Head_Pos_Time'; % name of figure to save
Walking = load([root '\Chirp_Walking_DATA.mat'],'PAT','HEAD','D','N', 'U'); % load data structure

FIG = figure (figNum); % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1100 800];
for kk = 1:Walking.N{1,1}
    for jj = 1:Walking.N{1,3}
        subplot(Walking.N{1,3},1,jj) ; hold on
        h.Trial = plot(Walking.HEAD.Time{kk}{jj},Walking.HEAD.Pos{kk}{jj},'Color',[0.5 0.5 0.5],'LineWidth',1); % trials
        for ii = 1:length(h.Trial)
            h.Trial(ii).Color(4) = 0.25;
        end
        h.Fly = plot(Walking.HEAD.FlyMed.Time{kk}(:,jj),Walking.HEAD.FlyMed.Pos{kk}(:,jj),'LineWidth',2); % flys
        h.Fly.Color(4) = 0.5;
    end
end

for jj = 1:Walking.N{1,3}
	subplot(Walking.N{1,3},1,jj) ; hold on
    title([num2str(Walking.U{1,3}{1}(jj)) '$^{\circ}$'],'Interpreter','latex','FontSize',15)
        h.patch = PlotPatch(Walking.HEAD.GrandMed.Pos(:,jj),Walking.HEAD.GrandSTD.Pos(:,jj),...
            Walking.HEAD.GrandMed.Time(:,jj),2,Walking.N{1,1},'k',[0.4 0.4 0.6],0.5,3); % all flys
    
     	xlim([0 20])
        ylim([-20 20])
        if jj==Walking.N{1,3}
            xlabel('Time (s)','Interpreter','latex','FontSize',15)
        end
        if jj~=Walking.N{1,3}
            xticks(0)
            xticklabels('')
        end
        ylabel('Head($^{\circ}$)','Interpreter','latex','FontSize',15)
end

saveas(FIG,[root '\' filename '.fig']); % save .fig file
print(gcf,[root '\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
disp('Saved to')
disp(root)
end
