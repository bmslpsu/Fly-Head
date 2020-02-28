function [] = MakeFig_SOS_CrossCorrelation()
%% MakeFig_SOS_CrossCorrelation:

root = 'H:\DATA\Rigid_Data\';

[FILES,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','on');
FILES = cellstr(FILES)';

nAmp = length(FILES);
Amp = nan(nAmp,1);
for ww = 1:nAmp
    filedata = textscan(FILES{ww}, '%s', 'delimiter', '_');
    Amp(ww) = str2double(filedata{1}{3});
end

HeadFree = cell(nAmp,1);
for ww = 1:nAmp
    HeadFree{ww} = load(fullfile(root,FILES{ww}),'TRIAL','GRAND','U','N');
end
%%
clearvars -except HeadFree Amp nAmp

s2h_idx = 5;
s2w_idx = 8;
h2w_idx = 7;

nFreq = HeadFree{1}.N{1,3};

lags = HeadFree{1}.GRAND{1,5}.Mean{2}{11}; % time lags
stim2head = cell(nFreq,nAmp);
stim2wing = cell(nFreq,nAmp);
head2wing = cell(nFreq,nAmp);
for ww = 1:nAmp
    nFly = HeadFree{ww}.N{1,1};
    for jj = 1:nFreq
        pp = 1;
        for kk = 1:nFly
            for ii = 1:size(HeadFree{ww}.TRIAL{kk,jj},1) % trials
                stim2head{jj,ww}(:,pp) = HeadFree{ww}.TRIAL{kk,jj}{ii,s2h_idx}.TimeDiff;
                stim2wing{jj,ww}(:,pp) = HeadFree{ww}.TRIAL{kk,jj}{ii,s2w_idx}.TimeDiff;
                head2wing{jj,ww}(:,pp) = HeadFree{ww}.TRIAL{kk,jj}{ii,h2w_idx}.TimeDiff;
                pp =  pp + 1;
            end
        end
    end
end

%% Head vs Wing
fig = figure (1); clf
set(fig,'Color','w','Units','inches', 'Position', [2 2 3 3], 'Name', 'Head vs Wing');
movegui(fig,'center')

ax = gobjects(nAmp,nFreq);

cc_idx = 5;
pp = 1;
for ww = 1:nAmp
    nFly = HeadFree{ww}.N{1,1};
    for jj = 1:nFreq
        ax(pp) = subplot(ceil(HeadFree{ww}.N{1,3}/1),1,pp) ; hold on
        for kk = 1:nFly
            for ii = 1:size(HeadFree{ww}.TRIAL{kk,jj},1) % trials
                CC      = HeadFree{ww}.TRIAL{kk,jj}{ii,cc_idx}.CrossCorr;
                MaxCC   = HeadFree{ww}.TRIAL{kk,jj}{ii,cc_idx}.MaxCC;
                TD      = HeadFree{ww}.TRIAL{kk,jj}{ii,cc_idx}.TimeDiff;
                
                plot(lags, CC, 'k', 'LineWidth', 0.5)
            	plot([TD TD], [-MaxCC MaxCC], 'r', 'LineWidth', 0.2)
            	plot(TD, MaxCC, '.', 'Color', 'r', 'MarkerSize', 10)
                %pause
             	%cla
            end
        end
        CC      = HeadFree{ww}.GRAND{jj,cc_idx}.Mean{1}{10};
        CC_STD  = HeadFree{ww}.GRAND{jj,cc_idx}.STD{1}{10};
        
        [~,max_idx] = max(abs(CC));
        TD = lags(max_idx);
        MaxCC = CC(max_idx);
        
        % MaxCC   = HeadFree{ww}.GRAND{jj,s2h_idx}.Mean{1}{12};
        % TD      = HeadFree{ww}.GRAND{jj,s2h_idx}.Mean{1}{13};
        
        PlotPatch(CC ,CC_STD, lags, 1, 1, 'c', [0.4 0.4 0.6], 0.5, 1);
        plot([TD TD TD], sort([MaxCC, 0 -MaxCC]), 'g', 'LineWidth', 1.5)
        plot(TD, MaxCC, '.', 'Color', 'g', 'MarkerSize', 10)
        pp = pp + 1;
    end
end
xlabel('Time Difference (s)')
ylabel('Cross Corelation')

set(ax, 'FontSize', 8, 'LineWidth', 1.5, 'XLim', 1*[-1 1])

end