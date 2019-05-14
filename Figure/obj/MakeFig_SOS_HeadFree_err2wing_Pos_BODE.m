function [FIG] = MakeFig_SOS_HeadFree_err2wing_Pos_BODE()
%% MakeFig_SOS_HeadFree_err2wing_Pos_BODE: BODE wing position for SOS
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

figNum = 1;
catIdx = 6;
xIdx = 1;

filename = 'SOS_HeadFree_err2wing_Pos_BODE'; % name of figure to save

FIG = figure (figNum); clf % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 800 600];
FIG.Name = filename;
movegui(FIG,'center')
hold on

% % Trials
% for kk = 1:HeadFree.N{1,1}
%     for ii = 1:size(HeadFree.TRIAL{kk},1)
%         subplot(2,1,1) ; hold on ; xlim([0.1 9]) ; ylim([0 2])
%             h.Trial = plot(HeadFree.TRIAL{kk}{ii,catIdx}.IOFreq,HeadFree.TRIAL{kk}{ii,catIdx}.IOBodeGain(:,xIdx),...
%                 '-','Color',[0.5 0.5 0.5 0.25],'LineWidth',1);
%         subplot(2,1,2) ; hold on ; xlim([0.1 9]) ; ylim(pi*[-1 1])
%             h.Trial = plot(HeadFree.TRIAL{kk}{ii,catIdx}.IOFreq,HeadFree.TRIAL{kk}{ii,catIdx}.IOBodePhaseDiff(:,xIdx),...
%                 '-','Color',[0.5 0.5 0.5 0.25],'LineWidth',1);       
%     end
% end

% Fly Stats
for kk = 1:HeadFree.N{1,1}
	subplot(2,1,1) ; hold on ; xlim([0.1 9]) ; ylim([0 2])
        h.Fly = plot(HeadFree.FLY{kk,catIdx}.Mean{4},HeadFree.FLY{kk,catIdx}.Mean{5}(:,xIdx),'-','LineWidth',2);
        h.Fly.Color(4) = 0.5;
   	subplot(2,1,2) ; hold on ; xlim([0.1 9]) ; ylim(pi*[-1 1])
        h.Fly = plot(HeadFree.FLY{kk,catIdx}.Mean{4},HeadFree.FLY{kk,catIdx}.CircMean{6}(:,xIdx),'-','LineWidth',2);
        h.Fly.Color(4) = 0.5;
end


% Grand Stats
subplot(2,1,1) ; hold on ; xlim([0 8.5]) ; ylim([0 3])
    errorbar(HeadFree.GRAND{1,catIdx}.Mean{2}{4},HeadFree.GRAND{1,catIdx}.Mean{2}{5}(:,xIdx),...
        2*HeadFree.GRAND{1,catIdx}.STD{2}{5}(:,xIdx),'-ok','LineWidth',3);
    ylabel('Head Gain ($^{\circ}/^{\circ}$)','Interpreter','latex','FontSize',15)
    xticks('')

subplot(2,1,2) ; hold on ; xlim([0 8.5]) ; ylim(pi*[-1 1])
    errorbar(HeadFree.GRAND{1,catIdx}.Mean{2}{4},HeadFree.GRAND{1,catIdx}.CircMean{9}{6}(:,xIdx),...
        2*HeadFree.GRAND{1,catIdx}.CircSTD{9}{6}(:,xIdx),'-ok','LineWidth',3);
%     plot([0 9],[0 0],'--g','LineWidth',2);
    ylabel('Head Phase (rad)','Interpreter','latex','FontSize',15)
    xlabel('Frequency (Hz)','Interpreter','latex','FontSize',15)

saveas(FIG,[root 'FIGURE\' filename '.fig']); % save .fig file
% print(gcf,[root 'FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
end