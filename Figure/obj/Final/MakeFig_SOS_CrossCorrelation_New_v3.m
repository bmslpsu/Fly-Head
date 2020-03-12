function [] = MakeFig_SOS_CrossCorrelation_New_v3()
%% MakeFig_SOS_CrossCorrelation_New_v3:

root = 'H:\DATA\Rigid_Data\';

[HeadFree,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head-free SOS data', root, 'MultiSelect','off');

HeadFree  = load(fullfile(root,HeadFree),'TRIAL','GRAND','U','N');

%%
clearvars -except HeadFree HeadFixed

Fs = 100;
scaleopt = 'normalized';
lags = HeadFree.GRAND{1,3}.Mean{2}{11}; % time lags
tthresh = 0.1;
trange = (lags >= -tthresh) & (lags <= tthresh);
lags = lags(trange);

ss = 1;
hh = 2;
ww = 2;

freqI = 1;

CC = cell(HeadFree.N.fly,3);
TD = cell(HeadFree.N.fly,3);
for kk = 1:HeadFree.N.fly % flys
    for ii = 1:size(HeadFree.TRIAL{kk,freqI},1) % trials
        stim = HeadFree.TRIAL{kk,freqI}{ii,ss}.X(:,1);
        head = HeadFree.TRIAL{kk,freqI}{ii,hh}.X(:,1);
        wing = HeadFree.TRIAL{kk,freqI}{ii,ww}.X(:,1);
        
        [cc,~,~,~] = CrossCorr(stim,head,Fs,scaleopt);
        CC{kk,1}(:,ii) = cc(trange);
        TD{kk,1}(:,ii) = get_time_diff(cc(trange),lags);
        
      	[cc,~,~,~] = CrossCorr(stim,wing,Fs,scaleopt);
        CC{kk,2}(:,ii) = cc(trange);
        TD{kk,2}(:,ii) = get_time_diff(cc(trange),lags);
        
      	[cc,~,~,~] = CrossCorr(head,wing,Fs,scaleopt);
        CC{kk,3}(:,ii) = cc(trange);
        TD{kk,3}(:,ii) = get_time_diff(cc(trange),lags);
        
%         for jj = 1:size(CC,2)
%             CC{kk,jj}(:,ii) = CC{kk,jj}(:,ii) / max(CC{kk,jj}(:,ii));
%         end
    end
end

TD = cellfun(@(x) median(x), TD, 'UniformOutput', true);
TD = median(TD,1);
TD = num2cell(TD);

Median_CC = cellfun(@(x) median(x,2), CC, 'UniformOutput', false);
% Median_TD = cellfun(@(x) median(x,1), TD, 'UniformOutput', true);
% Median_TD_All = median(Median_TD,1);
% STD_TD_All = std(Median_TD,[],1);

ngroup = size(Median_CC,2);
Median_CC_grouped = cell(1,ngroup);
for jj = 1:ngroup
    Median_CC_grouped{jj} = [Median_CC{:,jj}];
end

Median_CC_all = cellfun(@(x) median(x,2), Median_CC_grouped, 'UniformOutput', false);
STD_CC_all    = cellfun(@(x) std(x,[],2), Median_CC_grouped, 'UniformOutput', false);

fig = figure (1); clf
set(fig,'Color','w','Units','inches', 'Position', [2 2 3 3], 'Name', 'Head vs Wing');
movegui(fig,'center')
colors = {[0 0.7 0], [1 0 0], [0.3 0.1 0.7]};
ax = subplot(1,1,1) ; hold on
xlabel('Time Difference (s)')
ylabel('Normalized Cross-Correlation')
% TD = cell(1,ngroup);
MCC = cell(1,ngroup);
for jj = 1
    % h.fly = plot(lags, Median_CC_grouped{jj}, 'Color', [colors{jj}, 0.2], 'LineWidth', 0.5);
    [h.std(jj),h.med(jj)] = PlotPatch(Median_CC_all{jj}, STD_CC_all{jj}, lags, 1, 1, colors{jj}, ...
        0.5*colors{jj}, 0.15, 2);
    h.std(jj).EdgeColor = colors{jj};
    
    [~,MCC{jj}] = get_time_diff(Median_CC_all{jj},lags);
    
    h.td(jj) = plot([TD{jj}, TD{jj}] , [0 MCC{jj}], '.-', 'Color', colors{jj}, ...
                'LineWidth', 1 ,'MarkerSize', 15);
end
legend(h.med, 'stim2head', 'stim2wing', 'head2wing', 'box', 'off')

uistack(h.med,'top')
uistack(h.td,'top')
% xlim([-0.06 0.02])
set(ax, 'LineWidth', 1.5, 'FontSize', 8, 'Box', 'on', 'YLim', [0 1.05])

%%

fly = 1;
trial = 4;
idx_1 = 2;
idx_2 = 3;

x1 = HeadFree.TRIAL{fly}{trial,idx_1}.Time;
y1 = HeadFree.TRIAL{fly}{trial,idx_1}.X(:,1);

x2 = HeadFree.TRIAL{fly}{trial,idx_2}.Time;
y2 = HeadFree.TRIAL{fly}{trial,idx_2}.X(:,1);

close all
[xcf,lags,bounds,h] = crosscorr(y1,y2,'NumLags',1999);
% crosscorr(y1,y2,'NumLags',1999)
td = lags/100;

figure (1); clf
plot(td,xcf,'k')


scaleopt = 'normalized';
[R,lags] = xcorr(y1,y2,scaleopt);

figure (2); clf
plot(lags,R,'b')

% scaleopt = 'biased';
% [R,lags] = xcorr(y1,y2,scaleopt);
% 
% figure (1); clf
% subplot(2,1,1) ; cla ; hold on
% plot(x1,y1,'b')
% plot(x2,y2,'r')
% 
% subplot(2,1,2) ; cla ; hold on
% plot(lags,R,'k')
end

function [td,maxcc] = get_time_diff(cc,lags)
    [~,midx] = max((cc));
    maxcc = cc(midx);
    td = lags(midx);
end