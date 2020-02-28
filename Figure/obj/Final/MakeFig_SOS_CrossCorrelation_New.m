function [] = MakeFig_SOS_CrossCorrelation_New()
%% MakeFig_SOS_CrossCorrelation_New:

root = 'H:\DATA\Rigid_Data\';

[HeadFree,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head-free SOS data', root, 'MultiSelect','off');

[HeadFixed,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head-fixed SOS data', root, 'MultiSelect','off');

HeadFree  = load(fullfile(root,HeadFree),'TRIAL','GRAND','U','N');
HeadFixed = load(fullfile(root,HeadFixed),'TRIAL','GRAND','U','N');


%%
clearvars -except HeadFree HeadFixed

lags = HeadFree.GRAND{1,5}.Mean{2}{11}; % time lags
tthresh = 0.15;
trange = (lags >= -tthresh) & (lags <= tthresh);
lags = lags(trange);

s2h = 5;
h2w = 7;
s2w = 8;
s2w_fixed = 3;

CC = cell(HeadFree.N.fly,3);
for kk = 1:HeadFree.N.fly % flys
    for ii = 1:size(HeadFree.TRIAL{kk},1) % trials
        CC{kk,1}(:,ii) = HeadFree.TRIAL{kk}{ii,s2h}.CrossCorr(trange);
        CC{kk,2}(:,ii) = HeadFree.TRIAL{kk}{ii,s2w}.CrossCorr(trange);
        CC{kk,3}(:,ii) = HeadFree.TRIAL{kk}{ii,h2w}.CrossCorr(trange);       
        for jj = 1:size(CC,2)
            CC{kk,jj}(:,ii) = CC{kk,jj}(:,ii) / max(CC{kk,jj}(:,ii));
        end
    end
end

CC_fixed = cell(HeadFixed.N.fly,1);
for kk = 1:HeadFixed.N.fly % flys
    for ii = 1:size(HeadFixed.TRIAL{kk},1) % trials
        CC_fixed{kk,1}(:,ii) = HeadFixed.TRIAL{kk}{ii,s2w_fixed}.CrossCorr(trange);
        CC_fixed{kk}(:,ii) = CC_fixed{kk}(:,ii) / max(CC_fixed{kk}(:,ii));
    end
end

TD = cellfun(@(x) get_time_diff(x,lags), CC, 'UniformOutput', false);
TD_fixed = cellfun(@(x) get_time_diff(x,lags), CC_fixed, 'UniformOutput', false);

Median_CC = cellfun(@(x) median(x,2), CC, 'UniformOutput', false);
Median_TD = cellfun(@(x) median(x,1), TD, 'UniformOutput', true);
Median_TD_All = median(Median_TD,1);
% STD_TD_All = std(Median_TD,[],1);

Median_CC_fixed = cellfun(@(x) median(x,2), CC_fixed, 'UniformOutput', false);
Median_TD_fixed = cellfun(@(x) median(x,1), TD_fixed, 'UniformOutput', true);
Median_TD_All_fixed = median(Median_TD_fixed,1);
% STD_TD_All_fixed = std(Median_TD_fixed,[],1);

ngroup = size(Median_CC,2);
Median_CC_grouped = cell(1,ngroup);
for jj = 1:ngroup
    Median_CC_grouped{jj} = [Median_CC{:,jj}];
end
Median_CC_grouped_fixed = [Median_CC{:}];

Median_CC_all = cellfun(@(x) median(x,2), Median_CC_grouped, 'UniformOutput', false);
STD_CC_all    = cellfun(@(x) std(x,[],2), Median_CC_grouped, 'UniformOutput', false);

Median_CC_all_fixed = median(Median_CC_grouped_fixed, 2);
STD_CC_all_fixed    = std(Median_CC_grouped_fixed, [], 2);

fig = figure (1); clf
set(fig,'Color','w','Units','inches', 'Position', [2 2 3 3], 'Name', 'Head vs Wing');
movegui(fig,'center')
colors = {'b','r',[0.3 0.1 0.7] ,'g'};
ax = subplot(1,1,1) ; hold on
xlabel('Time Difference (s)')
ylabel('Normalized Cross-Correlation')
for jj = 1:ngroup
    h.fly = plot(lags, Median_CC_grouped{jj}, 'Color', [colors{jj}, 0.2], 'LineWidth', 0.5);
    [h.std(jj),h.med(jj)] = PlotPatch(Median_CC_all{jj}, STD_CC_all{jj}, lags, 0, 1, colors{jj}, ...
        [0.4 0.4 0.6], 0.2, 2);
    h.std(jj).EdgeColor = colors{jj};
    
    h.td(jj) = plot([Median_TD_All(jj), Median_TD_All(jj)] , [0 1], '.-', 'Color', colors{jj}, ...
                'LineWidth', 1 ,'MarkerSize', 15);
end

[h.std(ngroup+1),h.med(ngroup+1)] = PlotPatch(Median_CC_all_fixed, STD_CC_all_fixed, ...
                    lags, 1, 1, colors{ngroup+1}, [0.4 0.4 0.6], 0.2, 2);
                
h.td(ngroup+1) = plot([Median_TD_All_fixed, Median_TD_All_fixed] , [0 1], '.-', 'Color', colors{ngroup+1}, ...
                'LineWidth', 1 ,'MarkerSize', 15);

uistack(h.td,'top')
set(ax, 'LineWidth', 1.5, 'FontSize', 8, 'Box', 'on', 'YLim', [0 1.05])
end

function [td,maxcc] = get_time_diff(cc,lags)
    [~,midx] = max(abs(cc));
    maxcc = cc(midx);
    td = lags(midx);
end