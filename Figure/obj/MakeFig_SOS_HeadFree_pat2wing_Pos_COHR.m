function [FIG] = MakeFig_SOS_HeadFree_pat2wing_Pos_COHR()
%% MakeFig_SOS_HeadFree_pat2wing_Pos_COHR: head coherence plot for SOS
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'F:\DATA\Rigid_Data\';

% Select files
[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

HeadFree = load(fullfile(root,FILE{1}),'TRIAL','FLY','GRAND','U','N'); % load data structure
%%
figNum = 1;
catIdx = 8;
xIdx = 1;

filename = 'SOS_HeadFree_pat2wing_Pos_COHR'; % name of figure to save

FIG = figure (figNum); clf % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1100 800];
movegui(FIG,'center')
FIG.Name = filename;
hold on

% % Trials
% for kk = 1:HeadFree.N{1,1}
%     for ii = 1: size(HeadFree.TRIAL{kk},1)
%         h.Trial = plot(HeadFree.TRIAL{kk}{ii,catIdx}.IOFreq,HeadFree.TRIAL{kk}{ii,catIdx}.IOCoherence(:,xIdx),...
%             '-o','Color',[0.5 0.5 0.5 0.25],'LineWidth',1);
%     end
% end

% Fly Stats
for kk = 1:HeadFree.N{1,1}
    h.Fly = plot(HeadFree.FLY{kk,catIdx}.Mean{4},HeadFree.FLY{kk,catIdx}.Mean{9}(:,xIdx),'-o','LineWidth',2);
    h.Fly.Color(4) = 0.5;
end

% Grand Stats
h.patch = PlotPatch(HeadFree.GRAND{1,catIdx}.Mean{2}{7}(:,xIdx),HeadFree.GRAND{1,catIdx}.STD{2}{7}(:,xIdx),...
    HeadFree.GRAND{1,catIdx}.Mean{2}{8}(:,xIdx),2,HeadFree.N{1,1},'k',[0.4 0.4 0.6],0.5,3); % all flys

errorbar(HeadFree.GRAND{1,catIdx}.Mean{2}{4},HeadFree.GRAND{1,catIdx}.Mean{2}{9}(:,xIdx),...
    2*HeadFree.GRAND{1,catIdx}.STD{2}{9}(:,xIdx),'-ok','LineWidth',3);

xlim([0 9])
ylim([0 1])
xlabel('Frequency (Hz)','Interpreter','latex','FontSize',15)
ylabel('Pattern to Wing Coherence','Interpreter','latex','FontSize',15)

saveas(FIG,[root 'FIGURE\' filename '.fig']); % save .fig file
% print(gcf,[root 'FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
end