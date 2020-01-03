function [FIG] = MakeFig_Ramp_Head_Saccade_Int_Vel()
%% MakeFig_Ramp_Head_Saccade_Int_Vel:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

load(fullfile(root,FILE{1}),'SACCADE','INTERVAL','SACD','Stim','U','N','I','TRIAL','FLY','GRAND');

CC = repmat(hsv(N{1,3}/2),2,1);
% Vel = 3.75*U{1,3}{1};
Vel = U{1,3}{1};
clms = N.vel/2;

% clearvars -except SACCADE INTERVAL SACD Stim U N I TRIAL FLY GRAND CC Vel

%% Inter-Saccade Velocity %%
FIG = figure (1) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 (2)*clms 4];
FIG.Name = 'Inter-Saccade Head Velocity';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
clear h ax
ax = gobjects(N{1,3},1);
hp = gobjects(N{1,3},1);
hold on
pLim = 1;
for jj = 1:N{1,3}
	ax(jj) = subplot(ceil(N{1,3}/clms),clms,jj) ; hold on
    ax(jj).Title.String = [num2str(Vel(jj)) ' (' char(176) '/s)'];
    ax(jj).Title.Color = 'k';
    
    % Only plot until a percentage of intervals are still active
    remv = isnan(INTERVAL.Head.Position{jj}(1,:));
    INT = INTERVAL.Head.Position{jj}(:,~remv);
    tLim = sum(isnan(INT),2)./(size(INT,2)-1);
    span = 1:length(tLim(tLim<pLim));
    
%     plot(INTERVAL.Head.Time{jj},INTERVAL.Head.Velocity{jj},'Color', [0.5*CC{jj} , 0.2])
    
	[h.std(jj),h.med(jj)] = PlotPatch(INTERVAL.HeadStats.Velocity(jj).Median(span), INTERVAL.HeadStats.Velocity(jj).STD(span), ...
        INTERVAL.HeadStats.Time(jj).Median(span), 1, N.fly, CC(jj,:), [0.7 0.7 0.7], 0.4, 1.5);
    
    hp(jj) = plot([0 10], repmat(Vel(jj),1,2), '--', 'Color', CC(jj,:), 'LineWidth', 1.5);
end
set(ax,'FontSize',8,'Color','w','YColor','k','XColor','k','XLim',[0 1.5],'YLim',150*[-1 1])
XLabelHC = get(ax, 'XLabel');
set([XLabelHC{:}], 'String', 'Time (s)')
YLabelHC = get(ax([1,clms+1]), 'YLabel');
set([YLabelHC{:}], 'String', ['Head Velocity (' char(176) '/s)'])
linkaxes(ax,'xy')


%% Inter-Saccade Velocity ALL %%
FIG = figure (2) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 6 4];
FIG.Name = 'Inter-Saccade Head Velocity';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
clear h ax
ax = subplot(1,1,1) ; hold on
hp = gobjects(N{1,3},1);
pLim = 0.7;
for jj = 1:N{1,3}
    % Only plot until a percentage of intervals are still active
    remv = isnan(INTERVAL.Head.Position{jj}(1,:));
    INT = INTERVAL.Head.Position{jj}(:,~remv);
    tLim = sum(isnan(INT),2)./(size(INT,2)-1);
    span = 1:length(tLim(tLim<pLim));
    
    % plot(INTERVAL.Head.Time{jj},INTERVAL.Head.Velocity{jj},'Color', [0.5*CC(jj,:) , 0.2])
    
	[h.std(jj),h.med(jj)] = PlotPatch(INTERVAL.HeadStats.Velocity(jj).Median(span), INTERVAL.HeadStats.Velocity(jj).STD(span), ...
        INTERVAL.HeadStats.Time(jj).Median(span), 1, N.fly, CC(jj,:), [0.7 0.7 0.7], 0.4, 1.5);
    
    hp(jj) = plot([0 10], repmat(Vel(jj),1,2), '--', 'Color', 0.7*CC(jj,:), 'LineWidth', 1.5);
end
set(ax,'FontSize',8,'Color','w','YColor','k','XColor','k','XLim',[0 1.5],'YLim',150*[-1 1])
XLabelHC = get(ax, 'XLabel');
set(XLabelHC, 'String', 'Time (s)')
YLabelHC = get(ax, 'YLabel');
set(YLabelHC, 'String', ['Head Velocity (' char(176) '/s)'])
delete(h.std)

%% Inter-Saccade Velocity Error %%
FIG = figure (3) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 (2)*clms 4];
FIG.Name = 'Inter-Saccade Head Velocity Error';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
clear h ax
ax = gobjects(N{1,3},1);
hp = gobjects(N{1,3},1);
hold on
tt = 0:0.1:10;
pLim = 0.70;
for jj = 1:N{1,3}
	ax(jj) = subplot(ceil(N{1,3}/clms),clms,jj) ; hold on
    ax(jj).Title.String = [num2str(Vel(jj)) ' (' char(176) '/s)'];
    ax(jj).Title.Color = 'k';
    
    % Only plot until a percentage of intervals are still active
    remv = isnan(INTERVAL.Head.Position{jj}(1,:));
    INT = INTERVAL.Head.Position{jj}(:,~remv);
    tLim = sum(isnan(INT),2)./(size(INT,2)-1);
    span = 1:length(tLim(tLim<pLim));
    
    % plot(INTERVAL.Head.Time{jj},INTERVAL.Head.Velocity_Error{jj},'Color', [0.5*CC(jj,:) , 0.2])
    
	[h.std(jj),h.med(jj)] = PlotPatch(INTERVAL.HeadStats.Velocity_Error(jj).Median(span), INTERVAL.HeadStats.Velocity_Error(jj).STD(span), ...
        INTERVAL.HeadStats.Time(jj).Median(span), 1, N.fly, CC(jj,:), [0.7 0.7 0.7], 0.4, 2);
    
    hp(jj) = plot([0 10], repmat(Vel(jj),1,2), '--', 'Color', CC(jj,:), 'LineWidth', 1.5);
    
    plot(tt, 0.*tt, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
end
set(ax,'FontSize',12,'Color','w','YColor','k','XColor','k','XLim',[0 1.2],'YLim',150*[-1 1])
XLabelHC = get(ax(8:9), 'XLabel');
set([XLabelHC{:}], 'String', 'Time (s)')
YLabelHC = get(ax([1,clms+1]), 'YLabel');
set([YLabelHC{:}], 'String', ['Head Velocity Error (' char(176) '/s)'])
set([YLabelHC{:}], 'String', ['Retinal Slip (Velocity Error) (' char(176) '/s)'])
% delete(h.std)

%% Inter-Saccade Velocity Error ALL %%
FIG = figure (4) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 6 4];
FIG.Name = 'Inter-Saccade Head Velocity';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
clear h ax
ax = subplot(1,1,1) ; hold on
hp = gobjects(N{1,3},1);
pLim = 0.7;
for jj = 1:N{1,3}
    % Only plot until a percentage of intervals are still active
    remv = isnan(INTERVAL.Head.Position{jj}(1,:));
    INT = INTERVAL.Head.Position{jj}(:,~remv);
    tLim = sum(isnan(INT),2)./(size(INT,2)-1);
    span = 1:length(tLim(tLim<pLim));
    
    % plot(INTERVAL.Head.Time{jj},INTERVAL.Head.Velocity_Error{jj},'Color', [0.5*CC(jj,:) , 0.2])
    
	[h.std(jj),h.med(jj)] = PlotPatch(INTERVAL.HeadStats.Velocity_Error(jj).Median(span), INTERVAL.HeadStats.Velocity_Error(jj).STD(span), ...
        INTERVAL.HeadStats.Time(jj).Median(span), 1, N.fly, CC(jj,:), [0.7 0.7 0.7], 0.4, 1.5);
    
    hp(jj) = plot([0 10], repmat(Vel(jj),1,2), '--', 'Color', 0.7*CC(jj,:), 'LineWidth', 1.5);
end
plot([0 10], [0 0], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
set(ax,'FontSize',8,'Color','w','YColor','k','XColor','k','XLim',[0 1.5],'YLim',150*[-1 1])
XLabelHC = get(ax, 'XLabel');
set(XLabelHC, 'String', 'Time (s)')
YLabelHC = get(ax, 'YLabel');
set(YLabelHC, 'String', ['Head Velocity Error (' char(176) '/s)'])
uistack(h.med,'top')
leg = legend(h.med(1:clms),string(Vel(1:clms)));
leg.Box = 'off';
leg.Title.String = ['Speed (' char(176) '/s)'];
delete(h.std)

%% Inter-Saccade Velocity Error %%
FIG = figure (5) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 (2)*clms 4];
FIG.Name = 'Inter-Saccade Head Velocity Error';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
clear h ax
ax = gobjects(N{1,3},1);
hp = gobjects(N{1,3},1);
hold on
tt = 0:0.1:10;
pLim = 0.70;
for jj = 1:N{1,3}
	ax(jj) = subplot(ceil(N{1,3}/clms),clms,jj) ; hold on
    ax(jj).Title.String = [num2str(Vel(jj)) ' (' char(176) '/s)'];
    ax(jj).Title.Color = 'k';
    
    % Only plot until a percentage of intervals are still active
    remv = isnan(INTERVAL.Head.Position{jj}(1,:));
    INT = INTERVAL.Head.Position{jj}(:,~remv);
    tLim = sum(isnan(INT),2)./(size(INT,2)-1);
    span = 1:length(tLim(tLim<pLim));
    
    % plot(INTERVAL.Head.Time{jj},INTERVAL.Head.Velocity_Error{jj},'Color', [0.5*CC(jj,:) , 0.2])
    
	[h.std(jj),h.med(jj)] = PlotPatch(INTERVAL.HeadStats.Velocity_IntError(jj).Median(span), ...
        INTERVAL.HeadStats.Velocity_IntError(jj).STD(span), ...
        INTERVAL.HeadStats.Time(jj).Median(span), 1, N.fly, CC(jj,:), [0.7 0.7 0.7], 0.4, 2);
        
    plot(tt, 0.*tt, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
end
set(ax,'FontSize',12,'Color','w','YColor','k','XColor','k','XLim',[0 1.5],'YLim',30*[-1 1])
XLabelHC = get(ax(8:9), 'XLabel');
set([XLabelHC{:}], 'String', 'Time (s)')
YLabelHC = get(ax([1,clms+1]), 'YLabel');
set([YLabelHC{:}], 'String', ['Head Velocity Error (' char(176) '/s)'])
set([YLabelHC{:}], 'String', ['Retinal Slip (Velocity Error) (' char(176) '/s)'])

%% Inter-Saccade Velocity Integrated Error ALL %%
FIG = figure (6) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 6 4];
FIG.Name = 'Inter-Saccade Head Velocity';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
clear h ax
ax = subplot(1,1,1) ; hold on
hp = gobjects(N{1,3},1);
pLim = 0.7;
for jj = 1:N{1,3}
    % Only plot until a percentage of intervals are still active
    remv = isnan(INTERVAL.Head.Position{jj}(1,:));
    INT = INTERVAL.Head.Position{jj}(:,~remv);
    tLim = sum(isnan(INT),2)./(size(INT,2)-1);
    span = 1:length(tLim(tLim<pLim));
    
    % plot(INTERVAL.Head.Time{jj},INTERVAL.Head.Velocity_Error{jj},'Color', [0.5*CC(jj,:) , 0.2])
    
	[h.std(jj),h.med(jj)] = PlotPatch(INTERVAL.HeadStats.Velocity_IntError(jj).Median(span), ...
        INTERVAL.HeadStats.Velocity_IntError(jj).STD(span), ...
        INTERVAL.HeadStats.Time(jj).Median(span), 1, N.fly, CC(jj,:), [0.7 0.7 0.7], 0.4, 1.5);
    
%     hp(jj) = plot([0 10], repmat(Vel(jj),1,2), '--', 'Color', 0.7*CC(jj,:), 'LineWidth', 1.5);
end
plot([0 10], [0 0], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
set(ax,'FontSize',8,'Color','w','YColor','k','XColor','k','XLim',[0 1.5],'YLim',30*[-1 1])
XLabelHC = get(ax, 'XLabel');
set(XLabelHC, 'String', 'Time (s)')
YLabelHC = get(ax, 'YLabel');
set(YLabelHC, 'String', ['Head Velocity Error (' char(176) '/s)'])
uistack(h.med,'top')
leg = legend(h.med(1:clms),string(Vel(1:clms)));
leg.Box = 'off';
leg.Title.String = ['Speed (' char(176) '/s)'];
delete(h.std)

%% SS Gain %%
[b,a] = butter(2,2.5/(200/2),'low');
wind = 10;

figure (4) ; clf ; hold on
Gain = cell(N{1,3},1);
SS = cell(N{1,3},1);
TimeSS = cell(N{1,3},1);
pLim = 0.70;
for jj = 1:N{1,3}
  	% Only plot until a percentage of intervals are still active
    remv = isnan(INTERVAL.Head.Position{jj}(1,:));
    INT = INTERVAL.Head.Position{jj}(:,~remv);
    tLim = sum(isnan(INT),2)./(size(INT,2)-1);
    span = 1:length(tLim(tLim<pLim));
    for kk = 1:size(INTERVAL.Head.Velocity{jj},2)
        if ~isnan(INTERVAL.Head.Velocity{jj}(1,kk))
            time = INTERVAL.Head.Time{jj}(span,kk);
            vel = INTERVAL.Head.Velocity{jj}(span,kk);
            time = time(~isnan(time));
            vel = vel(~isnan(vel));
            svel = filtfilt(b,a,vel);
            if Vel(jj)>0
                [ss,peak_idx] = max(svel);
            elseif Vel(jj)<0
                [ss,peak_idx] = min(svel);
            end
        	peak_time = time(peak_idx);
            
            try
                % ss = mean(svel(peak_idx-wind:peak_idx+wind));
                ss = mean(svel);
            catch
                
            end
            
            if peak_time<10
                SS{jj}(end+1,1) = ss; 
                Gain{jj}(end+1,1) = ss/Vel(jj);
                TimeSS{jj}(end+1,1) = peak_time;
            end
                        
%             plot(time,vel)
%             plot(time,svel,'LineWidth',1.5)
%             plot(time,Vel(jj)*ones(length(time),1), '-k','LineWidth',1.5)
%             plot(peak_time,ss,'g.','MarkerSize',10)
%             plot([0 time(end)],[0 0], '--', 'Color', [0.5 0.5 0.5])
%             ylim(1.2*abs(Vel(jj))*[-1 1])
%             
%             pause
%             cla
        end
    end    
end

velIdx = num2cell((1:N.vel)');
G = cellfun(@(x,y) y*(x./x), Gain, velIdx, 'UniformOutput',  false);
G = cat(1,G{:});
G(G>clms) = G(G>clms) - clms;
Gain = cat(1,Gain{:});
SS = cat(1,SS{:});
TimeSS = cat(1,TimeSS{:});

DATA = {Gain,SS,TimeSS};
YY = {'Gain','Steady State (°/s)','SS Time (s)'};
% CC = prism(length(DATA));

FIG = figure (4) ; clf
rows = ceil(length(DATA)/3);
FIG.Units = 'inches';
FIG.Position = [2 2 2.5*3 2*rows];
FIG.Name = 'Inter-Saccade Head Velocity';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
clear h ax
ax = gobjects(length(DATA),1);
for jj = 1:length(DATA)
    ax(jj) = subplot(rows,3,jj); axis tight
    if jj==2
        data = abs(DATA{jj});
    else
        data = DATA{jj};
    end
    bx = boxplot(data, G,'Labels',{Vel(1:length(unique(G)))},'Width',0.5,'Symbol','','Whisker',2);
    xlabel(['Stimulus Speed (' char(176) '/s)'])
    ylabel(YY{jj})
    box off

    h = get(bx(5,:),{'XData','YData'});
    for kk = 1:size(h,1)
       patch(h{kk,1},h{kk,2},CC(kk,:));
    end

    set(findobj(ax(jj),'tag','Median'), 'Color', 'w','LineWidth',1.5);
    set(findobj(ax(jj),'tag','Box'), 'Color', 'none');
    set(findobj(ax(jj),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
    set(findobj(ax(jj),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
    ax(jj).Children = ax(jj).Children([end 1:end-1]);
end
set(ax,'FontSize',8)
ax(1).YLim = [-0.2 1.1];
ax(2).YLim = 80*[0 1];
ax(3).YLim = [0 1];

%% SS Gain MEAN %%
[b,a] = butter(2,2.5/(200/2),'low');
wind = 10;

figure (4) ; clf ; hold on
Gain = cell(N{1,3},1);
SS = cell(N{1,3},1);
TimeSS = cell(N{1,3},1);
pLim = 0.5;
for jj = 1:N{1,3}
	% Only plot until a percentage of intervals are still active
    remv = isnan(INTERVAL.Head.Position{jj}(1,:));
    INT = INTERVAL.Head.Position{jj}(:,~remv);
    tLim = sum(isnan(INT),2)./(size(INT,2)-1);
    span = 1:length(tLim(tLim<pLim));
    TIME = INTERVAL.HeadStats.Time(jj).Median(span);
    VEL_STD = INTERVAL.HeadStats.Velocity(jj).STD(span);
    VEL_MEAN = INTERVAL.HeadStats.Velocity(jj).Median(span);
    
    
    if Vel(jj)>0
        [ss,peak_idx] = max(VEL_MEAN);
    elseif Vel(jj)<0
        [ss,peak_idx] = min(VEL_MEAN);
    end
    peak_time = TIME(peak_idx);

    try
        % ss = mean(VEL_MEAN(peak_idx-wind:peak_idx+wind));
        ss = mean(VEL_MEAN);
    catch

    end

    if peak_time<10
        SS{jj}(end+1,1) = ss; 
        Gain{jj}(end+1,1) = ss/Vel(jj);
        TimeSS{jj}(end+1,1) = peak_time;
    end

%     % plot(TIME,vel)
%     plot(TIME,VEL_MEAN,'LineWidth',1.5)
%     plot(TIME,Vel(jj)*ones(length(TIME),1), '-k','LineWidth',1.5)
%     plot(peak_time,ss,'g.','MarkerSize',10)
%     plot([0 TIME(end)],[0 0], '--', 'Color', [0.5 0.5 0.5])
%     ylim(1.2*abs(Vel(jj))*[-1 1])
% 
%     pause
%     cla
end


%%
velIdx = num2cell((1:N.vel)');
G = cellfun(@(x,y) y*(x./x), Gain, velIdx, 'UniformOutput',  false);
G = cat(1,G{:});
G(G>clms) = G(G>clms) - clms;
Gain = cat(1,Gain{:});
SS = cat(1,SS{:});
TimeSS = cat(1,TimeSS{:});

DATA = {Gain,SS,TimeSS};
YY = {'Gain','Steady State (°/s)','SS Time (s)'};
CC = prism(length(DATA));

FIG = figure (4) ; clf
rows = ceil(length(DATA)/3);
FIG.Units = 'inches';
FIG.Position = [2 2 2.5*3 2*rows];
FIG.Name = 'Inter-Saccade Head Velocity';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
clear h ax
ax = gobjects(length(DATA),1);
for jj = 1:length(DATA)
    ax(jj) = subplot(rows,3,jj); axis tight
    if jj==2
        data = abs(DATA{jj});
    else
        data = DATA{jj};
    end
    bx = boxplot(data, G,'Labels',{Vel(1:length(unique(G)))},'Width',0.5,'Symbol','','Whisker',2);
    xlabel(['Stimulus Velocity (' char(176) '/s)'])
    ylabel(YY{jj})
    box off

    h = get(bx(5,:),{'XData','YData'});
    for k = 1:size(h,1)
       patch(h{k,1},h{k,2},CC(jj,:));
    end

    set(findobj(ax(jj),'tag','Median'), 'Color', 'w','LineWidth',1.5);
    set(findobj(ax(jj),'tag','Box'), 'Color', 'none');
    set(findobj(ax(jj),'tag','Upper Whisker'), 'Color', 'k');
    ax(jj).Children = ax(jj).Children([end 1:end-1]);
end
set(ax,'FontSize',8)
ax(1).YLim = [-1 2];
ax(2).YLim = 80*[0 1];
ax(3).YLim = [0 1];
end