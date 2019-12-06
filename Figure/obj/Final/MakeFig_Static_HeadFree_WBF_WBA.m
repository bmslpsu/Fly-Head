function [FIG] = MakeFig_Static_HeadFree_WBF_WBA()
%% MakeFig_Static_HeadFree_WBF_WBA:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

[Free,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head-free trials', root, 'MultiSelect','off');
Free = cellstr(Free)';

% [Fixed,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
%     'Select head-fixed trials', root, 'MultiSelect','off');
% Fixed = cellstr(Fixed)';

HeadFree  = load(fullfile(root,Free{1}),'TRIAL','GRAND','FLY','U','N');
% HeadFixed = load(fullfile(root,Fixed{1}),'TRIAL','GRAND','FLY','U','N');

%%
clearvars -except HeadFree HeadFixed
xIdx = 1;
head_free = 3;
head_fixed = 2;
catIdx = 3;

filename = 'Static_HeadFree_WBF_WBA';

nCond = HeadFree.N{1,3};

WBF = cell(nCond,1);

for jj = 1:nCond
    pp = 1;
    for kk = 1:HeadFree.N{1,1}
        for ii = 1:size(HeadFree.TRIAL{kk,jj},1)
           	WBF{jj,1}(pp,1) = kk;
            WBF{jj,1}(pp,2) = ii;
            WBF{jj,1}(pp,3) = false;
            WBF{jj,1}(pp,4) = median(HeadFree.TRIAL{kk,jj}{ii,head_free}.WBF(:,xIdx));
            pp = pp + 1;
        end
    end
end

bxdata = cat(1,WBF{:});
catg = nan(size(bxdata,1),1);
start = 1;
for jj = 1:nCond
 	nfly = size(WBF{jj},1);
    catg(start:(start+nfly-1)) = jj;
    start = length(catg(~isnan(catg))) + 1;
end


%% WBF Box Plot
FIG = figure (1); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 4 2];
FIG.Name = filename;
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')


% SOS_WBF = splitvars(array2table(DATA));
% SOS_WBF.Properties.VariableNames = {'Fly','Trial','Free_Fixed','WBF'};
% guide = 'Fixed is "true", Free is "False"';

label = string(HeadFree.U{1,3}{1});
label(1) = 'All Off';
label(5) = 'All On';
label(6) = 'Random';
bx = boxplot(bxdata(:,4),catg,'Labels',label,'Width',0.5,'Symbol','','Whisker',2);
% bx = boxplot(bxdata(:,4));

ax = gca;
ax.FontSize = 8;
ax.YLim = [200 300];
ylabel('WBF (Hz)','FontSize',8)
h = get(bx(5,:),{'XData','YData'});
CC = hsv(HeadFree.N{1,3});
for kk = 1:size(h,1)
   patch(h{kk,1},h{kk,2},CC(kk,:));
end
set(findobj(gcf,'tag','Median'), 'Color', 'w','LineWidth',1);
set(findobj(gcf,'tag','Box'), 'Color', 'k');
set(findobj(gcf,'tag','Upper Whisker'), 'Color', 'k','LineStyle','-','LineWidth',1);
set(findobj(gcf,'tag','Lower Whisker'), 'Color', 'k','LineStyle','-','LineWidth',1);
ax.Children = ax.Children([end 1:end-1]);
box off
%% WBF Raw
FIG = figure (1); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 6 3];
FIG.Name = filename;
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')

ax = gca; hold on
ax.FontSize = 12;
ax.XLabel.String = 'Time (s)';
ax.XLabel.FontSize = 14;
ax.YLabel.String = 'WBF (Hz)';
ax.YLabel.FontSize = 14;
ax.YLim = [180 280];

cList = lines(HeadFree.N{1,3});
clear h
pp = 1;
for jj = 1:HeadFree.N{1,3}
    for kk = 1:HeadFree.N{1,1}
        for ii = 1:size(HeadFree.TRIAL{kk,jj},1)
            h.trial = plot(HeadFree.TRIAL{kk,jj}{ii,catIdx}.Time, HeadFree.TRIAL{kk,jj}{ii,catIdx}.WBF,...
                'Color',cList(jj,:),'LineWidth',0.5);
            h.trial.Color(4) = 0.4;
            pp = pp + 1;
        end
        try
            h.fly(kk) = plot(HeadFree.FLY{jj}{kk,catIdx}.Mean{5}, HeadFree.FLY{jj}{kk,catIdx}.Mean{26},...
                'Color','r','LineWidth',1);
            catch
        end
    end
	h.grand(jj) = plot(HeadFree.GRAND{jj,catIdx}.Median{2}{5}, HeadFree.GRAND{jj,catIdx}.Median{2}{26},...
        'Color','k','LineWidth',2);
end
uistack(h.fly,'top')
uistack(h.grand,'top')

legend([h.trial(1),h.fly(1),h.grand(1)],'trials','fly means','grand mean')
%% WBA Raw
FIG = figure (2); clf
FIG.Color = 'w';
FIG.Position = [100 100 800 400];
FIG.Name = filename;
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')

ax = gca; hold on
ax.FontSize = 12;
ax.XLabel.String = 'Time (s)';
ax.XLabel.FontSize = 14;
ax.YLabel.String = 'WBA (Hz)';
ax.YLabel.FontSize = 14;
% ax.YLim = [200 280];

cList = prism(HeadFree.N.Amp);
clear h
for jj = 1:HeadFree.N.Amp
    for kk = 1:HeadFree.N.fly
        for ii = 1:size(HeadFree.TRIAL{kk,jj},1)
            h.trial = plot(HeadFree.TRIAL{kk,jj}{ii,catIdx}.Time, HeadFree.TRIAL{kk,jj}{ii,catIdx}.WBA(:,3),...
                'Color',cList(jj,:),'LineWidth',0.5);
            h.trial.Color(4) = 0.2;
        end
    end
	h.grand(jj) = plot(HeadFree.GRAND{jj,catIdx}.Mean{1}{5}, HeadFree.GRAND{jj,catIdx}.Mean{1}{27}(:,3),...
        'Color','k','LineWidth',3);
end
uistack(h.grand,'top')
end