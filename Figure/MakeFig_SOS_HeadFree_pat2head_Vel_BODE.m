function [] = MakeFig_SOS_HeadFree_pat2head_Vel_BODE()
%% MakeFig_SOS_HeadFree_pat2head_Vel_BODE: BODE head velocity for SOS
%   INPUTS:
%       -
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE INPUT %
root = 'F:\EXPERIMENTS\Experiment_SOS\DATA\';
figNum = 1;
%---------------------------------------------------------------------------------------------------------------------------------
filename = 'SOS_HeadFree_pat2head_Vel_BODE'; % name of figure to save
HeadFree = load([root 'SOS_HeadFree_DATA.mat'],'TRIAL','FLY','GRAND','U','N'); % load data structure

FIG = figure (figNum); % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1100 800];
FIG.Name = filename;
hold on
catIdx = 5;
xIdx = 2;

% Trials
for kk = 1:HeadFree.N{1,1}
    for ii = 1: size(HeadFree.TRIAL{kk},1)
        h.Trial = plot(HeadFree.TRIAL{kk}{ii,catIdx}.IOFreq,HeadFree.TRIAL{kk}{ii,catIdx}.IOBodeGain(:,xIdx),...
            '-o','Color',[0.5 0.5 0.5 0.25],'LineWidth',1);
    end
end

% Fly Stats
for kk = 1:HeadFree.N{1,1}
    h.Fly = plot(HeadFree.FLY{kk,catIdx}.Mean{4},HeadFree.FLY{kk,catIdx}.Mean{5}(:,xIdx),'-o','LineWidth',2);
    h.Fly.Color(4) = 0.5;
end

% Grand Stats
% h.patch = PlotPatch(HeadFree.GRAND{1,catIdx}.Mean{2}{5}(:,xIdx),HeadFree.GRAND{1,catIdx}.STD{2}{5}(:,xIdx),...
%     HeadFree.GRAND{1,catIdx}.Mean{2}{4},2,HeadFree.N{1,1},'k',[0.4 0.4 0.6],0.5,3); % all flys

errorbar(HeadFree.GRAND{1,catIdx}.Mean{2}{4},HeadFree.GRAND{1,catIdx}.Mean{2}{5}(:,xIdx),...
    2*HeadFree.GRAND{1,catIdx}.STD{2}{5}(:,xIdx),'-ok','LineWidth',3);

xlim([0 9])
ylim([0 1.2])
xlabel('Frequency (Hz)','Interpreter','latex','FontSize',15)
ylabel('Head Gain (${\circ}/{\circ}$)','Interpreter','latex','FontSize',15)

saveas(FIG,[root 'FIGURE\' filename '.fig']); % save .fig file
print(gcf,[root 'FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
disp('Saved to')
disp(root)
end