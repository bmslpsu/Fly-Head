function [] = MakeFig_SOS_Correlation_New()
%% MakeFig_SOS_Correlation:

root = 'H:\DATA\Rigid_Data\';

[HeadFree,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head-free SOS data', root, 'MultiSelect','off');

[HeadFixed,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head-free SOS data', root, 'MultiSelect','off');

HeadFree  = load(fullfile(root,HeadFree),'TRIAL','GRAND','U','N');
HeadFixed = load(fullfile(root,HeadFixed),'TRIAL','GRAND','U','N');


%%
clearvars -except HeadFree HeadFixed

stimIdx = 1;
headIdx = 2;
wingIdx = 3;
errIdx  = 4;
wingIdx_fixed = 2;
xIdx = 1;

STIM = cell(HeadFree.N.fly,1);
HEAD = cell(HeadFree.N.fly,1);
WING = cell(HeadFree.N.fly,1);
ERR  = cell(HeadFree.N.fly,1);

for kk = 1:HeadFree.N.fly % flys
    for ii = 1:size(HeadFree.TRIAL{kk},1) % trials
        STIM{kk}(:,ii) = HeadFree.TRIAL{kk}{ii,stimIdx}.X(:,xIdx);
        HEAD{kk}(:,ii) = HeadFree.TRIAL{kk}{ii,headIdx}.X(:,xIdx);
        WING{kk}(:,ii) = HeadFree.TRIAL{kk}{ii,wingIdx}.X(:,xIdx);
        ERR {kk}(:,ii) = HeadFree.TRIAL{kk}{ii,errIdx }.X(:,xIdx);
    end
end

STIM_fixed = cell(HeadFixed.N.fly,1);
WING_fixed = cell(HeadFixed.N.fly,1);
for kk = 1:HeadFixed.N.fly % flys
    for ii = 1:size(HeadFixed.TRIAL{kk},1) % trials        
        STIM_fixed{kk}(:,ii) = HeadFixed.TRIAL{kk}{ii,stimIdx}.X(:,xIdx);
        WING_fixed{kk}(:,ii) = HeadFixed.TRIAL{kk}{ii,wingIdx_fixed}.X(:,xIdx);
    end
end

stim2head = cell(HeadFree.N.fly,1);
stim2wing = cell(HeadFree.N.fly,1);
head2wing = cell(HeadFree.N.fly,1);
err2wing  = cell(HeadFree.N.fly,1);
for kk = 1:HeadFree.N.fly % flys
    for ii = 1:size(HeadFree.TRIAL{kk},1) % trials
        [stim2head{kk}(1,ii), stim2head{kk}(2,ii)] = corr(STIM{kk}(:,ii), HEAD{kk}(:,ii));
        [stim2wing{kk}(1,ii), stim2wing{kk}(2,ii)] = corr(STIM{kk}(:,ii), WING{kk}(:,ii));
        [head2wing{kk}(1,ii), head2wing{kk}(2,ii)] = corr(HEAD{kk}(:,ii), WING{kk}(:,ii));
        [err2wing{kk}(1,ii),  err2wing{kk}(2,ii) ] = corr(ERR{kk}(:,ii),  WING{kk}(:,ii));
        % [r,m,b] = regression(STIM{kk}(:,ii), HEAD{kk}(:,ii),'one');
    end
end

stim2wing_fixed = cell(HeadFixed.N.fly,1);
for kk = 1:HeadFixed.N.fly % flys
    for ii = 1:size(HeadFixed.TRIAL{kk},1) % trials
        [stim2wing_fixed{kk}(1,ii), stim2wing_fixed{kk}(2,ii)] = ...
            corr(STIM_fixed{kk}(:,ii), WING_fixed{kk}(:,ii));
    end
end

stim2head = cellfun(@(x) mean(x,2)', stim2head, 'UniformOutput', false);
stim2wing = cellfun(@(x) mean(x,2)', stim2wing, 'UniformOutput', false);
head2wing = cellfun(@(x) mean(x,2)', head2wing, 'UniformOutput', false);
err2wing  = cellfun(@(x) mean(x,2)', err2wing,  'UniformOutput', false);
stim2wing_fixed  = cellfun(@(x) mean(x,2)', stim2wing_fixed, 'UniformOutput', false);

stim2head = cat(1,stim2head{:});
stim2wing = cat(1,stim2wing{:});
head2wing = cat(1,head2wing{:});
err2wing  = cat(1,err2wing{:});
stim2wing_fixed = cat(1,stim2wing_fixed{:});

ALL = {stim2head, stim2wing, head2wing, err2wing, stim2wing_fixed};
MED = cellfun(@(x) median(x,1), ALL, 'UniformOutput', false);
MED = cat(1,MED{:});

names = {'stim2head', 'stim2wing', 'head2wing', 'err2wing', 'stim2wing_fixed'};
nameIdx = 1:length(ALL);
G = cellfun(@(x,y) y*ones(size(x,1),1), ALL, num2cell(nameIdx), 'UniformOutput', false);
G = cat(1,G{:});
ALL = cat(1,ALL{:});
colors = hsv(length(nameIdx));


FIG = figure (1) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 5 3];
movegui(FIG,'center')
ax(1) = subplot(1,1,1) ; hold on
bx = boxplot(ALL(:,1), G, 'Labels', names, 'Width', 0.5, 'Symbol', '+', 'Whisker', 2);
ylabel('Pearson Coefficent')
ylim([-1 1])
    h = get(bx(5,:),{'XData','YData'});
    for kk = 1:size(h,1)
       patch(h{kk,1},h{kk,2}, colors(kk,:));
    end
    plot([0,nameIdx,nameIdx(end)+1],zeros(2 +length(nameIdx),1), '--', ...
        'Color', [0.5 0.5 0.5], 'LineWidth', 0.75)

    set(findobj(ax(1),'tag','Median'), 'Color', 'w','LineWidth',1.5);
    set(findobj(ax(1),'tag','Box'), 'Color', 'none');
    set(findobj(ax(1),'tag','Upper Whisker'), 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
    set(findobj(ax(1),'tag','Lower Whisker'), 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
    ax(1).Children = ax(1).Children([end 1:end-1]);
    
set(ax, 'FontSize', 8,'LineWidth', 1.5, 'Box', 'off')
end