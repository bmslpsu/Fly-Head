function [FIG] = MakeFig_SOS_HeadFree_pat_Pos_Time()
%% MakeFig_SOS_HeadFree_pat_Pos_Time: time domain head position plot for SOS
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE INPUT %
root = 'H:\DATA\Rigid_Data\';

% Select files
[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

HeadFree = load(fullfile(root,FILE{1}),'TRIAL','FLY','GRAND','U','N'); % load data structure

figNum = 1;
catIdx = 1;
xIdx = 1;

filename = 'SOS_HeadFree_pat_Pos_Time'; % name of figure to save

FIG = figure (figNum); clf % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1100 800];
movegui(FIG,'center')
FIG.Name = filename;
hold on

% Trials
for kk = 1:1
    for ii = 1:1
        h.Trial = plot(HeadFree.TRIAL{kk}{ii,catIdx}.Time,HeadFree.TRIAL{kk}{ii,catIdx}.X(:,xIdx),...
            'Color',[0.5 0.5 0.5 1],'LineWidth',2);
    end
end

% % Fly Stats
% for kk = 1:HeadFree.N{1,1}
%     h.Fly = plot(HeadFree.FLY{kk,catIdx}.Mean{5},HeadFree.FLY{kk,catIdx}.Mean{6}(:,1),'LineWidth',2);
%     h.Fly.Color(4) = 0.5;
% end
% 
% % Grand Stats
% h.patch = PlotPatch(HeadFree.GRAND{1,catIdx}.Mean{2}{6}(:,1),HeadFree.GRAND{1,catIdx}.STD{2}{6}(:,1),...
%     HeadFree.GRAND{1,catIdx}.Mean{2}{5},2,HeadFree.N{1,1},'k',[0.4 0.4 0.6],0.5,3); % all flys

xlim([0 20])
ylim([-15 15])
xlabel('Time (s)','Interpreter','latex','FontSize',15)
ylabel('Pattern($^{\circ}$)','Interpreter','latex','FontSize',15)

saveas(FIG,[root 'FIGURE\' filename '.fig']); % save .fig file
%print(gcf,[root 'FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
end