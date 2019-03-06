function [] = MakeFig_Chirp_HeadFree_HeadvWing_Pos_Time(root,figNum)
%% MakeFig_Chirp_HeadFree_Head_Pos_Time:
%   INPUTS:
%       root: root directory containing data structure
%       figNum: figure #
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
% % EXAMPLE INPUT %
root = 'E:\EXPERIMENTS\Experiment_ChirpLog_HeadFree\DATA\';
figNum = 1;
%---------------------------------------------------------------------------------------------------------------------------------
filename = 'Chirp_HeadFree_HeadvWing_Pos_Time'; % name of figure to save
HeadFree = load([root 'Chirp_HeadFree_DATA.mat'],'PAT','WING','HEAD','n','unq'); % load data structure

FIG = figure (figNum); % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1400 600];
Fs = 200;
T = 20;
for jj = 1:HeadFree.n.Amp
    subplot(1,HeadFree.n.Amp,jj) ; hold on
    title([num2str(HeadFree.unq.Amp(jj)) '$^{\circ}$'],'Interpreter','latex','FontSize',15)
    xData = HeadFree.HEAD.ALL.Pos{jj}(1*Fs:T*Fs,:);
    yData = HeadFree.WING.ALL.Pos{jj}(1*Fs:T*Fs,:);
    [h] = scatplot(xData(:),yData(:));
    axis([-20 20 -8 8])
	delete(h.cb)
    xlabel('Head ($^{\circ}$)','Interpreter','latex','FontSize',15)
    if jj==1
        ylabel('Wing (V)','Interpreter','latex','FontSize',15)
    else
        yticks(0) ; yticklabels('')
    end

    % Calculate linear best fit
    [r,m,b] = regression(xData(:),yData(:),'one');
    text(20,7,['r =' num2str(r)]);
    xFit = linspace(-20,20,4000);
    yFit = m*xFit + b;
    plot(xFit,yFit,'r','LineWidth',5)
    
end

saveas(FIG,[root 'FIGURE\' filename '.fig']); % save .fig file
print(gcf,[root 'FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
disp('Saved to')
disp(root)
end
