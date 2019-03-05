function [] = MakeFig_ALL_HeadFixed_Wing_Pos_BODE()
%% MakeFig_ALL_HeadFixed_Wing_Pos_BODE:
%   INPUTS:
%       root: root directory containing data structure
%       figNum: figure #
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
root        = 'E:\DATA\';
filename    = 'ALL_HeadFixed_Wing_Pos_BODE'; % name of figure to save
HeadFree    = load([root 'Chirp_HeadFree_DATA.mat'],'PAT','WING','HEAD','n','unq');
HeadFixed   = load([root 'Chirp_HeadFixed_DATA.mat'],'PAT','WING','n','unq');
Linear      = load([root 'ChirpLinear15_HeadFree_DATA.mat'],'PAT','WING','n','unq');
Replay      = load([root 'ChirpReplay15_HeadFixed_DATA.mat'],'PAT','WING','n','unq');

FIG = figure (1); % figure handle
FIG.Color = 'w';
set(gcf,'Position',[100 100 1000 700])
lineColor = {'b','r','c','g'};
faceColor = [0.5 0.5 0.5];
alphaVal = 0.5;
lineWidth = 2;
AMP = 3;
legLabel = {'Head Free','Head Fixed','Replay','Linear'};
% MAGNITUDE
subplot(2,1,1) ; hold on ; title([num2str(HeadFixed.unq.Amp(3)) '$^{\circ}$'],'Interpreter','latex','FontSize',15)
	ylabel('Wing Gain($V/{\circ}$)','Interpreter','latex','FontSize',15)
    h.HeadFree  = PlotPatch(HeadFree.WING.GrandMean.GAIN(:,AMP),HeadFree.WING.GrandSTD.GAIN(:,AMP),...
                    HeadFree.WING.GrandMean.Freq(:,AMP),2,HeadFree.n.Fly,lineColor{1},faceColor,alphaVal,lineWidth);
    h.HeadFixed = PlotPatch(HeadFixed.WING.GrandMean.GAIN(:,AMP),HeadFixed.WING.GrandSTD.GAIN(:,AMP),...
                    HeadFixed.WING.GrandMean.Freq(:,AMP),2,HeadFixed.n.Fly,lineColor{2},faceColor,alphaVal,lineWidth);
    h.Replay    = PlotPatch(Replay.WING.GrandMean.GAIN(:,1),Replay.WING.GrandSTD.GAIN(:,1),...
                    Replay.WING.GrandMean.Freq(:,1),2,Replay.n.Fly,lineColor{3},faceColor,alphaVal,lineWidth);
    h.Linear    = PlotPatch(Linear.WING.GrandMean.GAIN(:,1),Linear.WING.GrandSTD.GAIN(:,1),...
                    Linear.WING.GrandMean.Freq(:,1),2,Linear.n.Fly,lineColor{4},faceColor,alphaVal,lineWidth);                
	xlim([0.1 12])
    ylim([0 0.35])
    xticks(0)
    xticklabels('')
    legend([h.HeadFree h.HeadFixed h.Replay h.Linear],legLabel)
% PHASE    
subplot(2,1,2) ; hold on
    xlabel('Frequency (Hz)','Interpreter','latex','FontSize',15)
    ylabel('Wing Phase(rad)','Interpreter','latex','FontSize',15)
    h.HeadFree  = PlotPatch(HeadFree.WING.GrandMean.PHASE(:,AMP),HeadFree.WING.GrandSTD.PHASE(:,AMP),...
                    HeadFree.WING.GrandMean.Freq(:,AMP),2,HeadFree.n.Fly,lineColor{1},faceColor,alphaVal,lineWidth);
    h.HeadFixed = PlotPatch(HeadFixed.WING.GrandMean.PHASE(:,AMP),HeadFixed.WING.GrandSTD.PHASE(:,AMP),...
                    HeadFixed.WING.GrandMean.Freq(:,AMP),2,HeadFixed.n.Fly,lineColor{2},faceColor,alphaVal,lineWidth);
    h.Replay    = PlotPatch(Replay.WING.GrandMean.PHASE(:,1),Replay.WING.GrandSTD.PHASE(:,1),...
                    Replay.WING.GrandMean.Freq(:,1),2,Replay.n.Fly,lineColor{3},faceColor,alphaVal,lineWidth);
    h.Linear    = PlotPatch(Linear.WING.GrandMean.PHASE(:,1),Linear.WING.GrandSTD.PHASE(:,1),...
                    Linear.WING.GrandMean.Freq(:,1),2,Linear.n.Fly,lineColor{4},faceColor,alphaVal,lineWidth);                
	xlim([0.1 12])
    xticks([0.1 2:2:12])
    ylim([-3 3])

%---------------------------------------------------------------------------------------------------------------------------------
saveas(FIG,[root 'FIGURE\' filename '.fig']); % save .fig file
print(gcf,[root 'FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
disp('Saved to')
disp(root)
end
