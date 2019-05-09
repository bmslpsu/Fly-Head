function [] = MakeFig_SOS_HeadFree_Head_Pos_Time()
%% MakeFig_SOS_HeadFree_Pos_Time: time domain head position plot for SOS
%   INPUTS:
%       -
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE INPUT %
root = 'F:\EXPERIMENTS\Experiment_SOS\DATA\';
figNum = 1;
%---------------------------------------------------------------------------------------------------------------------------------
filename = 'SOS_HeadFree_Head_Pos_Time'; % name of figure to save
HeadFree = load([root 'SOS_HeadFree_DATA.mat'],'TRIAL','FLY','GRAND','U','N'); % load data structure

FIG = figure (figNum); % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1100 800];
hold on
catIdx = 2;

% Trials
for kk = 1:HeadFree.N{1,1}
    for ii = size(HeadFree.TRIAL{kk},1)
        h.Trial = plot(HeadFree.TRIAL{kk}{ii,catIdx}.Time,HeadFree.TRIAL{kk}{ii,catIdx}.X(:,1),...
            'Color',[0.5 0.5 0.5 0.25],'LineWidth',1);
    end
end

% Fly Stats
for kk = 1:HeadFree.N{1,1}
    h.Fly = plot(HeadFree.FLY{kk,catIdx}.Mean{5},HeadFree.FLY{kk,catIdx}.Mean{6}(:,1),'LineWidth',2);
    h.Fly.Color(4) = 0.5;
end

% Grand Stats
h.patch = PlotPatch(HeadFree.GRAND{1,catIdx}.Mean{2}{6}(:,1),HeadFree.GRAND{1,catIdx}.STD{2}{6}(:,1),...
    HeadFree.GRAND{1,catIdx}.Mean{2}{5},2,HeadFree.N{1,1},'k',[0.4 0.4 0.6],0.5,3); % all flys

xlim([0 20])
ylim([-15 15])
xlabel('Time (s)','Interpreter','latex','FontSize',15)
ylabel('Head($^{\circ}$)','Interpreter','latex','FontSize',15)

saveas(FIG,[root 'FIGURE\' filename '.fig']); % save .fig file
print(gcf,[root 'FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
disp('Saved to')
disp(root)
end