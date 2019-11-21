function [FIG] = MakeFig_Ramp_Head_Saccade_Int_Pos_v2()
%% MakeFig_Ramp_Head_Saccade_Int_Pos_v2:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

load(fullfile(root,FILE{1}),'INTERVAL','Stim','U','N','T');

clearvars -except SACCADE INTERVAL SACD Stim U N I TRIAL FLY GRAND CC Vel

clms = N.vel/2;
CC = repmat(hsv(clms),2,1);
Vel = U.vel{1};

%% Inter-Saccade Time %%
interval_times = cell(N.vel,1);
for jj = 1:N.vel
	TIME = INTERVAL.Head.Time{jj};
    POS = INTERVAL.Head.Position{jj};
    n_int = size(TIME,2);
    for kk = 1:n_int % don't use first interval
        valid_times = TIME(~isnan(POS(:,kk)),kk);
        if isempty(valid_times)
            interval_times{jj}(end+1,1) = nan;
        else
            interval_times{jj}(end+1,1) = valid_times(end);
        end
    end
end

med_times = cellfun(@(x) nanmedian(x), interval_times, 'UniformOutput',  true);
std_times = cellfun(@(x) nanstd(x), interval_times, 'UniformOutput',  true);

velIdx = num2cell((1:N.vel)');
G = cellfun(@(x,y) y*(x./x), interval_times, velIdx, 'UniformOutput',  false);
G = cat(1,G{:});
G(G>clms) = G(G>clms) - clms;
interval_times_all = cat(1,interval_times{:});

interval_cut.pos = cell(N{1,3},1);
lim = med_times + 1*std_times;
for jj = 1:N{1,3}
    for kk = 1:size(INTERVAL.Head.Position{jj},2)
        if interval_times{jj}(kk)>lim(jj)
            interval_cut.pos{jj}(:,kk) = INTERVAL.Head.Position{jj}(:,kk);
        else
            interval_cut.pos{jj}(:,kk) = nan*INTERVAL.Head.Position{jj}(:,kk);
        end
    end
end
interval_cut_median.pos = cellfun(@(x) nanmean(x,2), interval_cut.pos, 'UniformOutput', false);
% nanIdx = cellfun(@(x) ~isnan(x), interval_cut_median.pos, 'UniformOutput', false);
% interval_cut_median.pos = cellfun(@(x,y) x(y), interval_cut_median.pos, nanIdx, 'UniformOutput', false);

%% Inter-Saccade Times %%
FIG = figure (1) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 (2/3)*clms 2];
FIG.Name = 'Inter-Saccade Times';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')

ax = gca;
bx = boxplot(interval_times_all,G,'Labels',{Vel(1:clms)},'Width',0.5,'Symbol','','Whisker',2);
xlabel(['Stimulus Velocity (' char(176) '/s)'])
ylabel('Interval Time (s)')
box off

h = get(bx(5,:),{'XData','YData'});
for kk = 1:size(h,1)
   patch(h{kk,1},h{kk,2},CC(kk,:));
end

set(findobj(ax,'tag','Median'), 'Color', 'w','LineWidth',1.5);
set(findobj(ax,'tag','Box'), 'Color', 'none');
set(findobj(ax,'tag','Upper Whisker'), 'Color', 'k');
ax.Children = ax.Children([end 1:end-1]);
ax.YLim = [0 3];

%%
FIG = figure (11) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 5 5];
FIG.Name = 'Inter-Saccade Times';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
ax = gobjects(N.vel,1);
shape = reshape(1:N.vel,N.vel/2,2);
for jj = 1:N.vel/2
    % shape(jj)
    ax(jj) = subplot(1,1,1) ; hold on
    histogram(interval_times{jj},100,'FaceColor',CC(jj,:))
end


%% Inter-Saccade Position %%
FIG = figure (2) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 6 4];
FIG.Name = 'Inter-Saccade Head Position';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
clear ax h

hp = gobjects(N{1,3},1);
ax = axes;
hold on

% pLim = 0.50;
Ts = 1/200;
for jj = 1:N{1,3}
    % Only plot until a percentage of intervals are still active
%     remv = isnan(INTERVAL.Head.Time{jj}(1,:));
%     INT = INTERVAL.Head.Time{jj}(:,~remv);
%     tLim = sum(isnan(INT),2)./(size(INT,2)-1);
%     span = 1:length(tLim(tLim<pLim));
    
    tlim = med_times(jj);
    stdlim = med_times(jj) + 1*std_times(jj);
    span_med = 1:round(tlim/Ts);
    span_std = 1:round(stdlim/Ts);
    
    plot(INTERVAL.Head.Time{jj}(:,:),INTERVAL.Head.Position{jj}(:,:),'Color', [0.7*CC(jj,:) , 0.2])
    
	[h.std(jj),h.med(jj)] = PlotPatch(INTERVAL.HeadStats.Position(jj).Median(span_med), INTERVAL.HeadStats.Position(jj).STD(span_med), ...
        INTERVAL.HeadStats.Time(jj).Median(span_med), 1, 1, CC(jj,:), [0.7 0.7 0.7], 0.4, 3);
    delete(h.std(jj))
    
 	[h.std(jj),h.m(jj)] = PlotPatch(INTERVAL.HeadStats.Position(jj).Median(span_std), INTERVAL.HeadStats.Position(jj).STD(span_std), ...
        INTERVAL.HeadStats.Time(jj).Median(span_std), 1, 1, CC(jj,:), [0.7 0.7 0.7], 0.4, 3);
    delete(h.m(jj))
    
	hp(jj) = plot(INTERVAL.HeadStats.Time(jj).Median(span_med), Stim(span_med,jj), '--', 'Color', 0.5*CC(jj,:), 'LineWidth', 1.5);
    
end
plot([0 10],[0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',1.5)
set(ax,'FontSize',8,'Color','w','YColor','k','XColor','k','XLim',[0 1.5])
ax.YLim = max(abs(ax.YLim))*[-1 1];
ax.YLim = 30*[-1 1];
ax.XLabel.String = 'Time (s)';
ax.YLabel.String = ['Head Angle (' char(176) ')'];
uistack(h.std,'top')
uistack(h.med,'top')
uistack(hp,'top')


leg = legend(h.med(1:clms),string(Vel(1:clms)));
leg.Box = 'off';
leg.Title.String = ['Speed (' char(176) '/s)'];

%% Inter-Saccade Position Error %%
FIG = figure (3) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 6 4];
FIG.Name = 'Inter-Saccade Head Position Error';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
clear ax h

ax = gca;
hp = gobjects(N{1,3},1);
hold on

pLim = 0.50;
for jj = 1:N{1,3}
    % Only plot until a percentage of intervals are still active
    remv = isnan(INTERVAL.Head.Time{jj}(1,:));
    INT = INTERVAL.Head.Time{jj}(:,~remv);
    tLim = sum(isnan(INT),2)./(size(INT,2)-1);
    span_med = 1:length(tLim(tLim<pLim));
    
    h.trial = plot(INTERVAL.Head.Time{jj},INTERVAL.Head.Position_Error{jj},'Color', [0.5*CC(jj,:),0.2]);
  	[h.std(jj),h.med(jj)] = PlotPatch(INTERVAL.HeadStats.Position_Error(jj).Median(span_med), INTERVAL.HeadStats.Position_Error(jj).STD(span_med), ...
        INTERVAL.HeadStats.Time(jj).Median(span_med), 1, 1, CC(jj,:), [0.7 0.7 0.7], 0.4, 3);
    
    hp(jj) = plot(INTERVAL.HeadStats.Time(jj).Median(span_med), Stim(span_med,jj), '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 2);
    
    plot([0 10],[0 0], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
    
%     plot([0 ; repmat(INTERVAL.HeadStats.Time(jj).Median(span(end)),2,1)],...
%         [0 , INTERVAL.HeadStats.Position_Error(jj).Median(span(end)) , 0],...
%         '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
    
% 	plot([repmat(INTERVAL.HeadStats.Time(jj).Median(span(end)),2,1)],...
%         [INTERVAL.HeadStats.Position_Error(jj).Median(span(end)) , 0],...
%         '--', 'Color', CC(jj,:), 'LineWidth', 2);
    
end
set(ax,'FontSize',8,'Color','w','YColor','k','XColor','k','XLim',[0 1.5])
ax.YLim = max(abs(ax.YLim))*[-1 1];
ax.YLim = 60*[-1 1];
ax.XLabel.String = 'Time (s)';
ax.YLabel.String = ['Head Position Error (' char(176) ')'];
uistack(h.std,'top')
uistack(h.med,'top')
uistack(hp,'top')

%% Inter-Saccade Integrated Position Error %%
FIG = figure (5) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 6 4];
FIG.Name = 'Inter-Saccade Head Position Error';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
clear ax h

ax = gca;
hold on
int_level = nan(N.vel,1);
t_level = nan(N.vel,1);
pLim = 0.5;
for jj = 1:N{1,3}
    % Only plot until a percentage of intervals are still active
    remv = isnan(INTERVAL.Head.Time{jj}(1,:));
    INT = INTERVAL.Head.Time{jj}(:,~remv);
    tLim = sum(isnan(INT),2)./(size(INT,2)-1);
    span_med = 1:length(tLim(tLim<pLim));
    
    h.trial = plot(INTERVAL.Head.Time{jj},INTERVAL.Head.Position_IntError{jj},'Color', [0.5*CC(jj,:),0.2]);
    
    int_level(jj) = INTERVAL.HeadStats.Position_IntError(jj).Median(span_med(end));
    t_level(jj) = INTERVAL.HeadStats.Time(jj).Median(span_med(end));
    
  	[h.std(jj),h.med(jj)] = PlotPatch(INTERVAL.HeadStats.Position_IntError(jj).Median(span_med), INTERVAL.HeadStats.Position_IntError(jj).STD(span_med), ...
        INTERVAL.HeadStats.Time(jj).Median(span_med), 1, 1, CC(jj,:), [0.7 0.7 0.7], 0.4, 3);
    
    h.int_level(jj) = plot([0 t_level(jj)], [int_level(jj) int_level(jj)], '--', 'Color', CC(jj,:), 'LineWidth', 1.5);
    h.int_marker(jj) = plot(t_level(jj), int_level(jj), '.k', 'MarkerSize', 20);
    
end
plot([0 10],[0 0], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
set(ax,'FontSize',8,'Color','w','YColor','k','XColor','k','XLim',[0 1.5])
ax.YLim = max(abs(ax.YLim))*[-1 1];
ax.YLim = 15*[-1 1];
ax.XLabel.String = 'Time (s)';
ax.YLabel.String = ['Head Integrated Position Error (' char(176) '*s)'];
uistack(h.std,'top')
uistack(h.med,'top')
uistack(h.int_level,'top')
uistack(h.int_marker,'top')

%% Inter-Saccade Position Error AREA %%
FIG = figure (4) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 6 4];
FIG.Name = 'Inter-Saccade Head Position Error';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
clear ax h

ax = gca;
hp = gobjects(N{1,3},1);
hold on

pLim = 0.5;
for jj = 1:N{1,3}
    % Only plot until a percentage of intervals are still active
    remv = isnan(INTERVAL.Head.Time{jj}(1,:));
    INT = INTERVAL.Head.Time{jj}(:,~remv);
    tLim = sum(isnan(INT),2)./(size(INT,2)-1);
    span_med = 1:length(tLim(tLim<pLim));
    
%     h.trial = plot(INTERVAL.Head.Time{jj},INTERVAL.Head.Position_Error{jj},'Color', [0.5*CC(jj,:),0.2]);
%   	[h.std(jj),h.med(jj)] = PlotPatch(INTERVAL.HeadStats.Position_Error(jj).Median(span), INTERVAL.HeadStats.Position_Error(jj).STD(span), ...
%         INTERVAL.HeadStats.Time(jj).Median(span), 1, 1, CC(jj,:), [0.7 0.7 0.7], 0.4, 3);
    
%     hp(jj) = plot(INTERVAL.HeadStats.Time(jj).Median(span), Stim(span,jj), '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 2);
    
    plot([0 10],[0 0], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
    
    plot([0 ; repmat(INTERVAL.HeadStats.Time(jj).Median(span_med(end)),2,1)],...
        [0 , INTERVAL.HeadStats.Position_Error(jj).Median(span_med(end)) , 0],...
        '-', 'Color', CC(jj,:), 'LineWidth', 2);
    
	plot([repmat(INTERVAL.HeadStats.Time(jj).Median(span_med(end)),2,1)],...
        [INTERVAL.HeadStats.Position_Error(jj).Median(span_med(end)) , 0],...
        '--', 'Color', CC(jj,:), 'LineWidth', 2);
    
end
set(ax,'FontSize',8,'Color','w','YColor','k','XColor','k','XLim',[0 1.5])
ax.YLim = max(abs(ax.YLim))*[-1 1];
ax.YLim = 60*[-1 1];
ax.XLabel.String = 'Time (s)';
ax.YLabel.String = ['Head Position Error (' char(176) ')'];
% uistack(h.std,'top')
% uistack(h.med,'top')
% uistack(hp,'top')

%% Inter-Saccade Position CUT %%
FIG = figure (10) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 6 4];
FIG.Name = 'Inter-Saccade Head Position';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
clear ax h

hp = gobjects(N{1,3},1);
ax = axes;
hold on

pLim = 1;
for jj = 1:N{1,3}
    % Only plot until a percentage of intervals are still active
    remv = isnan(INTERVAL.Head.Time{jj}(1,:));
    INT = INTERVAL.Head.Time{jj}(:,~remv);
    tLim = sum(isnan(INT),2)./(size(INT,2)-1);
    span_med = 1:length(tLim(tLim<pLim));
    
    plot(INTERVAL.Head.Time{jj},interval_cut.pos{jj},'Color', [0.7*CC(jj,:) , 0.2])

    plot(INTERVAL.HeadStats.Time(jj).Median,interval_cut_median.pos{jj},'Color',CC(jj,:),'LineWidth',3)
    
%  	[h.std(jj),h.med(jj)] = PlotPatch(INTERVAL.HeadStats.Position(jj).Median(span), INTERVAL.HeadStats.Position(jj).STD(span), ...
%         INTERVAL.HeadStats.Time(jj).Median(span), 1, 1, CC(jj,:), [0.7 0.7 0.7], 0.4, 3);
    
	hp(jj) = plot(INTERVAL.HeadStats.Time(jj).Median(span_med), Stim(span_med,jj), '--', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5);
    
end
plot([0 10],[0 0],'--','Color',[0.5 0.5 0.5],'LineWidth',1.5)
set(ax,'FontSize',8,'Color','w','YColor','k','XColor','k','XLim',[0 1.5])
ax.YLim = max(abs(ax.YLim))*[-1 1];
ax.YLim = 30*[-1 1];
ax.XLabel.String = 'Time (s)';
ax.YLabel.String = ['Head Angle (' char(176) ')'];
% uistack(h.std,'top')
% uistack(h.med,'top')
uistack(hp,'top')


% leg = legend(h.med(1:clms),string(Vel(1:clms)));
% leg.Box = 'off';
% leg.Title.String = ['Speed (' char(176) '/s)'];







end