function [FIG] = MakeFig_SOS_HeadFree_head_wing_TimeDiff()
%% MakeFig_SOS_HeadFree_head_wing_TimeDiff: 
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%

root = 'H:\DATA\Rigid_Data\';

% Select files (Head fixed: SOS_HeadFixed_DATA_06-12-2019)
[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

HeadFree = load(fullfile(root,FILE{1}),'TRIAL','GRAND','FLY','U','N');

%%
headIdx = 3;
wingIdx = 3;

xIdx = 1;

filename = 'SOS_HeadFree_pat2head_BODE';

FIG = figure (1); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 4 2];
FIG.Name = filename;
movegui(FIG,'center')
hold on

ax1 = subplot(1,1,1) ; hold on
    ax1.FontSize = 8;
    ax1.XLim = [0 10];
    ax1.YLim = [-100 100];
    ax1.XLabel.FontSize = 8;
    ax1.XLabel.Color = 'k';
    ax1.XLabel.String = 'Frequency (Hz)';
 	ax1.YLabel.String = 'Time Difference (ms)';
    ax1.YLabel.FontSize = ax1.XLabel.FontSize;
    
%     HFreq           = HeadFree.GRAND{1,headIdx}.Mean{2}{4};
%     HPhase          = rad2deg(HeadFree.GRAND{1,headIdx}.CircMean{7}{6}(:,xIdx));
%     HPhase_STD      = rad2deg(HeadFree.GRAND{1,headIdx}.CircSTD{7}{6}(:,xIdx));
%     HTimeDiff       = 1000*(HPhase/360).*(1./HFreq);
%     HTimeDiff_STD   = 1000*(HPhase_STD/360).*(1./HFreq);
    
%     [~,h.gain] = PlotPatch(HTimeDiff, HTimeDiff_STD, HFreq, 1, HeadFree.N{1,1}, 'b', [0.4 0.4 0.6], 0.5, 2);
% 	h.gain.Marker = '.';
%     h.gain.MarkerSize = 20;
    
	WFreq           = HeadFree.GRAND{1,wingIdx}.Mean{2}{4};
    WPhase          = rad2deg(HeadFree.GRAND{1,wingIdx}.CircMean{7}{6}(:,xIdx));
    cnd             = find(WPhase>30);
    WPhase(cnd)   	= WPhase(cnd) -360;
    
    WPhase_STD      = rad2deg(HeadFree.GRAND{1,wingIdx}.CircSTD{7}{6}(:,xIdx));
    WTimeDiff       = 1000*(WPhase/360).*(1./WFreq);
    WTimeDiff_STD   = 1000*(WPhase_STD/360).*(1./WFreq);
    
    [~,h.phase] = PlotPatch(WTimeDiff, WTimeDiff_STD, WFreq, 1, HeadFree.N{1,1}, 'r', [0.4 0.4 0.6], 0.5, 2);
	h.phase.Marker = '.';
    h.phase.MarkerSize = 20;
    
%% Box Plot acorss frequencies

wfreq = HeadFree.GRAND{1,wingIdx}.Mean{2}{4};
wphs  = nan(length(wfreq),HeadFree.N.fly);
wtd   = nan(length(wfreq),HeadFree.N.fly);
for kk = 1:HeadFree.N.fly
	wphase          = rad2deg(HeadFree.FLY{kk,wingIdx}.CircMean{6}(:,xIdx));
    
    cnd             = wphase > -130 & any(((1:length(wphase))'==[4,5]),2);
    wphase(cnd)   	= wphase(cnd) - 360;
%     
%     cnd             = (wphase < -300) & ((1:length(wphase))'==5);
%     wphase(cnd)   	= wphase(cnd) + 360;
    
    
    wphs(:,kk)      = wphase;
 	wtd(:,kk)     	= 1000*(wphase/360).*(1./wfreq);
end

% cla ; hold on
% plot(WFreq,wphs)
% plot(WFreq, rad2deg(circ_mean(deg2rad(wphs),[],2)),'k','LineWidth',2)
% [~,h.phase] = PlotPatch(WPhase, WPhase_STD, WFreq, 1, HeadFree.N{1,1}, 'r', [0.4 0.4 0.6], 0.5, 2);
% h.phase.Marker = '.';
% h.phase.MarkerSize = 20;

FIG = figure (2) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 2 2];
movegui(FIG,'center')
clear ax
ww = 1;
ax(ww) = subplot(1,1,ww); axis tight

bx = boxplot(wtd(:), 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
ylabel('Time Difference (ms)')

h = get(bx(5,:),{'XData','YData'});
for kk = 1:size(h,1)
   patch(h{kk,1},h{kk,2},'k');
end

set(findobj(ax(ww),'tag','Median'), 'Color', 'w','LineWidth',1.5);
set(findobj(ax(ww),'tag','Box'), 'Color', 'none');
set(findobj(ax(ww),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
set(findobj(ax(ww),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
ax(ww).Children = ax(ww).Children([end 1:end-1]);
ylim([-150 0])
set(ax, 'LineWidth', 1, 'Box', 'on')


end