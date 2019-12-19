function [FIG] = MakeFig_SOS_Walking_pat2head_Bode_all()
%% MakeFig_SOS_HeadFree_pat_head_Time:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'E:\Walking_Experiments';

[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

Walking = load(fullfile(root,FILE{1}),'TRIAL','GRAND', 'FLY','U','N', 'T');
%%
figNum = 1;
catIdx = 4;
xIdx = 1;

filename = 'SOS_Walking_all_head_Time';

FIG = figure (figNum); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 8 2];
movegui(FIG,'center')
FIG.Name = filename;

ax.L = gca; hold on
ax.L.YColor = 'k';
ax.L.FontSize = 12;
ax.L.YLabel.String = ['Magnitude'];
ax.L.YLabel.FontSize = 12;
% ax.L.YLim =[0 10];
% ax.L.YTick = 15*[-1 0 1];
ax.L.XLabel.String = 'Frequency';
ax.L.XLabel.Color = 'k';
ax.L.XLabel.FontSize = 12;
ax.L.XLim = [0 12];
ax.L.XTick = [1, 3.1, 5.3, 7.4, 9.6];

% Trials
% for kk = 1:Walking.N{1,1}
%     for jj = 1:Walking.T{kk,2}
%         subplot(2,1,1) ; hold on ; xlim([0.1 9]) ; ylim([0 2])
%         plot(Walking.TRIAL{kk}{jj,catIdx}.IOFreq,Walking.TRIAL{kk}{jj,catIdx}.IOBodeGain(:,xIdx),'Color',[0.5 0.5 0.5 0.25],'LineWidth',1);
%         subplot(2,1,2) ; hold on ; xlim([0.1 9]) ; ylim(pi*[-1 1])
%         tphase = rad2deg(Walking.TRIAL{kk}{jj,catIdx}.IOBodePhaseDiff(:,xIdx));
%         tphase(4) = tphase(4)-360;
%         plot(Walking.TRIAL{kk}{jj,catIdx}.IOFreq,tphase,'Color',[0.5 0.5 0.5 0.25],'LineWidth',1);      
%     end
% end

% Fly Stats
% for kk = 1:Walking.N{1,1}
% 	subplot(2,1,1) ; hold on ; xlim([0.1 9]) ; ylim([0 2])
%     plot(Walking.FLY{kk,catIdx}.Mean{4},Walking.FLY{kk,catIdx}.Mean{5}(:,xIdx),'LineWidth',1);
%     subplot(2,1,2) ; hold on ; xlim([0.1 9]) ; ylim(pi*[-1 1])
%     fphase = rad2deg(Walking.FLY{kk,catIdx}.Mean{6}(:,xIdx));
%     fphase(4) = fphase(4)-360;
%     plot(Walking.FLY{kk,catIdx}.Mean{4},fphase,'LineWidth',1);
% end

% Grand Stats
ax1 = subplot(2,1,1) ; hold on
    ax1.FontSize = 12;
    ax1.XLim = [0 9.6]; 
    ax1.XTick = [1, 3.1, 5.3, 7.4, 9.6];
    ax1.YLim = [0 1];
    ax1.XTickLabel = '';
	ax1.XLabel.String = 'Frequency (Hz)';
    ax1.XLabel.FontSize = 14;
    ax1.XLabel.Color = 'w';
 	ax1.YLabel.String = ['Head Gain (' char(176) '/' char(176) ')'];
    ax1.YLabel.FontSize = ax1.XLabel.FontSize;
    
    h.patch = PlotPatch(Walking.GRAND{1,catIdx}.Mean{2}{5}(:,xIdx),Walking.GRAND{1,catIdx}.STD{2}{5}(:,xIdx),...
        Walking.GRAND{1,catIdx}.Mean{2}{4},2,Walking.N{1,1},'k',[0.4 0.4 0.6],0.5,3);
    
ax2 = subplot(2,1,2) ; hold on
    ax2.FontSize = ax1.FontSize;
    ax2.XLim = ax1.XLim;
    ax2.XTick = ax1.XTick;
    ax2.YLim = [-300 60];
   	ax2.XLabel.String = 'Frequency (Hz)';
    ax2.XLabel.FontSize = 14;
    ax2.XLabel.Color = 'k';
 	ax2.YLabel.String = ['Head Phase (' char(176) ')'];
    ax2.YLabel.FontSize = ax1.XLabel.FontSize;
    %ax2.YTick = -180:30:180;
    gphase = rad2deg(Walking.GRAND{1,catIdx}.CircMean{7}{6}(:,xIdx));
    gphase(4) = gphase(4)-360;
    gphase(5) = gphase(5)-360;
    gstdphase = rad2deg(Walking.GRAND{1,catIdx}.CircSTD{7}{6}(:,xIdx));
    h.patch = PlotPatch(gphase,gstdphase,Walking.GRAND{1,catIdx}.Mean{2}{4},2,Walking.N{1,1},'k',[0.4 0.4 0.6],0.5,3);
    plot(ax1.XLim,[0 0],'--', 'Color', [0.4660 0.6740 0.1880],'LineWidth',2);


end