function [FIG] = MakeFig_SOS_HeadFree_Variance()
%% MakeFig_SOS_HeadFree_Variance: 
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

HeadFree = load(fullfile(root,FILE{1}),'TRIAL','FLY','GRAND','U','N');
%%
headIdx = 5;
wingIdx = 6;

xIdx = 1;

filename = 'SOS_HeadFree_pat2head_BODE';

FIG = figure (1); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = 1*[1 1 4 4];
FIG.Name = filename;
movegui(FIG,'center')
CC = jet(HeadFree.N.fly);

ax(1) = subplot(2,1,1); hold on
ax(1).YLim = [0 1];
ax(1).XLim = [0 10];
clear h
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

%%
HEAD = nan(HeadFree.N.fly,2);
for kk = 1:HeadFree.N.fly
	HEAD(kk,1) = mean(HeadFree.FLY{kk,headIdx}.STD{5}(:,xIdx));
    HEAD(kk,2) = mean(rad2deg(HeadFree.FLY{kk,headIdx}.CircSTD{6}(:,xIdx)));
end


end