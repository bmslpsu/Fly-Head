function [FIG] = MakeFig_SOS_HeadFree_err2wing_BODE()
%% MakeFig_SOS_HeadFree_err2wing_BODE: BODE head position for SOS
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

% Select files
[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

HeadFree = load(fullfile(root,FILE{1}),'TRIAL','FLY','GRAND','U','N'); % load data structure
%%
figNum = 1;
catIdx = 6;
xIdx = 1;

filename = 'SOS_HeadFree_err2wing_BODE'; % name of figure to save

FIG = figure (figNum); clf % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 600 600];
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

% % Fly Stats
% for kk = 1:HeadFree.N{1,1}
% 	subplot(2,1,1) ; hold on ; xlim([0.1 9]) ; ylim([0 2])
%         h.Fly = plot(HeadFree.FLY{kk,catIdx}.Mean{4},HeadFree.FLY{kk,catIdx}.Mean{5}(:,xIdx),'-','LineWidth',2);
%         h.Fly.Color(4) = 0.5;
%    	subplot(2,1,2) ; hold on ; xlim([0.1 9]) ; ylim(pi*[-1 1])
%         h.Fly = plot(HeadFree.FLY{kk,catIdx}.Mean{4},HeadFree.FLY{kk,catIdx}.CircMean{6}(:,xIdx),'-','LineWidth',2);
%         h.Fly.Color(4) = 0.5;
% end

% Grand Stats
ax1 = subplot(2,1,1) ; hold on
    ax1.FontSize = 12;
    ax1.XLim = [0 10]; 
    ax1.YLim = [0 0.3];
    ax1.XTickLabel = '';
% 	ax1.XLabel.String = 'Frequency (Hz)';
    ax1.XLabel.FontSize = 14;
    ax1.XLabel.Color = 'w';
 	ax1.YLabel.String = ['Error to Wing Gain (V/' char(176) ')'];
    ax1.YLabel.FontSize = ax1.XLabel.FontSize;
    
    errorbar(HeadFree.GRAND{1,catIdx}.Mean{2}{4},HeadFree.GRAND{1,catIdx}.Mean{2}{5}(:,xIdx),...
        2*HeadFree.GRAND{1,catIdx}.STD{2}{5}(:,xIdx),'-or','LineWidth',3);
    
ax2 = subplot(2,1,2) ; hold on
    ax2.FontSize = ax1.FontSize;
    ax2.XLim = ax1.XLim;
    ax2.YLim = [-300 60];
   	ax2.XLabel.String = 'Frequency (Hz)';
    ax2.XLabel.FontSize = 14;
    ax2.XLabel.Color = 'k';
 	ax2.YLabel.String = ['Error to Wing Phase (' char(176) ')'];
    ax2.YLabel.FontSize = ax1.XLabel.FontSize;
    ax2.YTick = -300:60:180;
    
    phase = HeadFree.GRAND{1,catIdx}.CircMean{9}{6}(:,xIdx);
    phase(phase>1) = phase(phase>1) - 2*pi;
    
    errorbar(HeadFree.GRAND{1,catIdx}.Mean{2}{4},rad2deg(phase),...
        rad2deg(2*HeadFree.GRAND{1,catIdx}.CircSTD{9}{6}(:,xIdx)),'-or','LineWidth',3);
    
    plot(ax1.XLim,[0 0],'--g','LineWidth',2);

end