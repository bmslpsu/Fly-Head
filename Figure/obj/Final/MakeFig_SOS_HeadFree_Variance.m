function [FIG] = MakeFig_SOS_HeadFree_Variance()
%% MakeFig_SOS_HeadFree_Variance: 
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';
% root = 'Q:\';
% Select files
[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

HeadFree = load(fullfile(root,FILE{1}),'TRIAL','FLY','GRAND','U','N');
%%
clearvars -except HeadFree

headIdx = 5;
wingIdx = 6;
xIdx = 1;

%% Head FRF
FIG = figure (1); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = 1.5*[1 1 4 4];
movegui(FIG,'center')
CC = jet(HeadFree.N.fly);
clear h ax

ax(1) = subplot(2,1,1); hold on
ax(1).YLim = [0 1];
ax(1).XLim = [0 10];
for kk = 1:HeadFree.N.fly
    for ii = 1:size(HeadFree.TRIAL{kk},1)
        plot(HeadFree.TRIAL{kk}{ii,headIdx}.IOFreq, HeadFree.TRIAL{kk}{ii,headIdx}.IOBodeGain(:,xIdx),...
            'Color',[0.5 0.5 0.5 0.3],'Marker','.','MarkerSize',5)
    end   
    [~,h.fly(kk)] = PlotPatch(HeadFree.FLY{kk,headIdx}.Mean{5}(:,xIdx), HeadFree.FLY{kk,headIdx}.STD{5}(:,xIdx), ...
        HeadFree.FLY{kk,headIdx}.Mean{4}, 1, HeadFree.N.fly, CC(kk,:), [0.4 0.4 0.6], 0.4, 1);
end
set(h.fly,'Marker','.','MarkerSize',15)
uistack(h.fly,'top')

[~,h.grand] = PlotPatch(HeadFree.GRAND{headIdx}.Mean{2}{5}(:,xIdx), HeadFree.GRAND{headIdx}.STD{2}{5}(:,xIdx), ...
    HeadFree.GRAND{headIdx}.Mean{2}{4}, 1, HeadFree.N.fly, 'k', [0.4 0.4 0.6], 0.5, 3);
set(h.grand,'Marker','.','MarkerSize',30)
ylabel(['Gain (' char(176) '/' char(176) ')'])


ax(2) = subplot(2,1,2); hold on
ax(2).YLim = [-120 60];
ax(2).XLim = [0 10];
clear h
for kk = 1:HeadFree.N.fly
    for ii = 1:size(HeadFree.TRIAL{kk},1)
        plot(HeadFree.TRIAL{kk}{ii,headIdx}.IOFreq, rad2deg(HeadFree.TRIAL{kk}{ii,headIdx}.IOBodePhaseDiff(:,xIdx)),...
            'Color',[0.5 0.5 0.5 0.3],'Marker','.','MarkerSize',5)
    end   
    [~,h.fly(kk)] = PlotPatch(rad2deg(HeadFree.FLY{kk,headIdx}.CircMean{6}(:,xIdx)), rad2deg(HeadFree.FLY{kk,headIdx}.CircSTD{6}(:,xIdx)), ...
        HeadFree.FLY{kk,headIdx}.Mean{4}, 1, HeadFree.N.fly, CC(kk,:), [0.4 0.4 0.6], 0.4, 1);
end
set(h.fly,'Marker','.','MarkerSize',15)
uistack(h.fly,'top')

[~,h.grand] = PlotPatch(rad2deg(HeadFree.GRAND{headIdx}.CircMean{2}{6}(:,xIdx)), rad2deg(HeadFree.GRAND{headIdx}.CircSTD{8}{6}(:,xIdx)), ...
    HeadFree.GRAND{headIdx}.Mean{2}{4}, 1, HeadFree.N.fly, 'k', [0.4 0.4 0.6], 0.5, 3);
set(h.grand,'Marker','.','MarkerSize',30)

xlabel('Frequency (Hz)')
ylabel(['Phase (' char(176) ')'])

set(ax,'FontSize',8)
linkaxes(ax,'x')

%% Wing FRF
FIG = figure (2); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = 1.5*[1 1 4 4];
movegui(FIG,'center')
CC = jet(HeadFree.N.fly);

clear h ax
ax(1) = subplot(2,1,1); hold on
ax(1).YLim = [0 0.4];
ax(1).XLim = [0 10];
for kk = 1:HeadFree.N.fly
    for ii = 1:size(HeadFree.TRIAL{kk},1)
        plot(HeadFree.TRIAL{kk}{ii,wingIdx}.IOFreq, HeadFree.TRIAL{kk}{ii,wingIdx}.IOBodeGain(:,xIdx),...
            'Color',[0.5 0.5 0.5 0.3],'Marker','.','MarkerSize',5)
    end   
    [~,h.fly(kk)] = PlotPatch(HeadFree.FLY{kk,wingIdx}.Mean{5}(:,xIdx), HeadFree.FLY{kk,wingIdx}.STD{5}(:,xIdx), ...
        HeadFree.FLY{kk,wingIdx}.Mean{4}, 1, HeadFree.N.fly, CC(kk,:), [0.4 0.4 0.6], 0.4, 1);
end
set(h.fly,'Marker','.','MarkerSize',15)
uistack(h.fly,'top')

[~,h.grand] = PlotPatch(HeadFree.GRAND{wingIdx}.Mean{2}{5}(:,xIdx), HeadFree.GRAND{wingIdx}.STD{2}{5}(:,xIdx), ...
    HeadFree.GRAND{wingIdx}.Mean{2}{4}, 1, HeadFree.N.fly, 'k', [0.4 0.4 0.6], 0.5, 3);
set(h.grand,'Marker','.','MarkerSize',30)
ylabel(['Gain (V/' char(176) ')'])


ax(2) = subplot(2,1,2); hold on
ax(2).YLim = [-360 120];
ax(2).XLim = [0 10];
clear h
for kk = 1:HeadFree.N.fly
    for ii = 1:size(HeadFree.TRIAL{kk},1)
        plot(HeadFree.TRIAL{kk}{ii,wingIdx}.IOFreq, rad2deg(HeadFree.TRIAL{kk}{ii,wingIdx}.IOBodePhaseDiff(:,xIdx)),...
            'Color',[0.5 0.5 0.5 0.3],'Marker','.','MarkerSize',5)
    end   
    [~,h.fly(kk)] = PlotPatch(rad2deg(HeadFree.FLY{kk,wingIdx}.CircMean{6}(:,xIdx)), rad2deg(HeadFree.FLY{kk,wingIdx}.CircSTD{6}(:,xIdx)), ...
        HeadFree.FLY{kk,wingIdx}.Mean{4}, 1, HeadFree.N.fly, CC(kk,:), [0.4 0.4 0.6], 0.4, 1);
end
set(h.fly,'Marker','.','MarkerSize',15)
uistack(h.fly,'top')

[~,h.grand] = PlotPatch(rad2deg(HeadFree.GRAND{wingIdx}.CircMean{2}{6}(:,xIdx)), rad2deg(HeadFree.GRAND{wingIdx}.CircSTD{8}{6}(:,xIdx)), ...
    HeadFree.GRAND{wingIdx}.Mean{2}{4}, 1, HeadFree.N.fly, 'k', [0.4 0.4 0.6], 0.5, 3);
set(h.grand,'Marker','.','MarkerSize',30)

xlabel('Frequency (Hz)')
ylabel(['Phase (' char(176) ')'])

set(ax,'FontSize',8)
linkaxes(ax,'x')

%% Head & Wings Variance
Head = nan(5,HeadFree.N.fly,2);
Wing = nan(5,HeadFree.N.fly,2);
for kk = 1:HeadFree.N.fly
	Head(:,kk,1) = HeadFree.FLY{kk,headIdx}.STD{5}(:,xIdx);
    Head(:,kk,2) = rad2deg(HeadFree.FLY{kk,headIdx}.CircSTD{6}(:,xIdx));
  	Wing(:,kk,1) = HeadFree.FLY{kk,wingIdx}.STD{5}(:,xIdx);
    Wing(:,kk,2) = rad2deg(HeadFree.FLY{kk,wingIdx}.CircSTD{6}(:,xIdx));
end
HEAD = mean(Head,1);
WING = mean(Wing,1);

FIG = figure (3); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = 1.5*[1 1 4 4];
movegui(FIG,'center')
bins = 8;
clear ax

ax(1) = subplot(2,2,1); hold on ; title('Head STD') ; xlabel(['Gain (' char(176) '/' char(176) ')'])
histogram(HEAD(:,:,1),bins,'FaceColor','b')

ax(2) = subplot(2,2,2); hold on ; title('Wing STD') ; xlabel(['Gain (V/' char(176) ')'])
histogram(WING(:,:,1),bins,'FaceColor','r')

ax(3) = subplot(2,2,3); hold on ; xlabel(['Phase (' char(176) ')'])
histogram(HEAD(:,:,2),bins,'FaceColor','b')

ax(4) = subplot(2,2,4); hold on ; xlabel(['Phase (' char(176) ')'])
histogram(WING(:,:,2),bins,'FaceColor','r')

set(ax,'FontSize',8,'YLim',[0 5],'Box','on')

%% Head & Wings Variance Combined
Head = nan(5,HeadFree.N.fly,2);
Wing = nan(5,HeadFree.N.fly,2);
for kk = 1:HeadFree.N.fly
	Head(:,kk,1) = HeadFree.FLY{kk,headIdx}.STD{5}(:,xIdx);
    Head(:,kk,2) = rad2deg(HeadFree.FLY{kk,headIdx}.CircSTD{6}(:,xIdx));
  	Wing(:,kk,1) = HeadFree.FLY{kk,wingIdx}.STD{5}(:,xIdx);
    Wing(:,kk,2) = rad2deg(HeadFree.FLY{kk,wingIdx}.CircSTD{6}(:,xIdx));
end
HEAD = mean(Head,1);
WING = mean(Wing,1);

FIG = figure (4); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 4 2];
movegui(FIG,'center')
edges = 5:5:45;
clear ax

ax(1) = subplot(1,1,1); hold on ; title('Phase STD') ; xlabel(['Phase (' char(176) ')']) ; ylabel('Count')
histogram(HEAD(:,:,2),edges,'FaceColor','b','FaceAlpha',1)
histogram(WING(:,:,2),edges,'FaceColor','r','FaceAlpha',1)
legend('Head','Wings','box','off')

set(ax,'FontSize',8,'YLim',[0 10],'Box','on')
end