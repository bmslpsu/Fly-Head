function [FIG] = MakeFig_Ramp_Head_Saccade()
%% MakeFig_Ramp_Head_Saccade:
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

CC = repmat({[1 0 0],[0 1 0],[0 0 1]},1,2);

Vel = 3.75*U{1,3}{1};
clearvars -except SACCADE INTERVAL SACD Stim U N I TRIAL FLY GRAND CC Vel
%% Saccade Position %%
FIG = figure (1) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 6 4];
FIG.Name = 'Saccade Position';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = axes;
for jj = 1:N{1,3}
    ax(jj) = subplot(ceil(N{1,3}/3),3,jj) ; hold on
    ax(jj).Title.String = [num2str(Vel(jj)) ' (' char(176) '/s)'];
    ax(jj).Title.Color = 'k';
    
    cent = max(SACCADE.cIdx{jj},[],'all');
    span = (cent-5):(cent+5);
    
    plot(1000*SACCADE.Head.Time{jj}, SACCADE.Head.Position{jj}, 'Color', 0.7*CC{jj})
%     plot(1000*SACCADE.HeadStats.Time(jj).Median(span),SACCADE.HeadStats.Position(jj).Median(span),'w','LineWidth',3)
    
    PlotPatch(SACCADE.HeadStats.Position(jj).Median(span), SACCADE.HeadStats.Position(1,1).STD(span), ...
        1000*SACCADE.HeadStats.Time(jj).Median(span), 3, N.fly, CC{jj}, [0.7 0.7 0.7], 0.4, 3);
    
    if jj<=N{1,3}/2
        ax(jj).XTickLabels = '';
    else
        ax(jj).XLabel.String = 'Time (ms)';
    end
    
    if (jj-1)/3 ~= round((jj-1)/3)
        ax(jj).YTickLabels = '';
    else
        ax(jj).YLabel.String = ['(' char(176) ')'];
    end
end
set(ax,'FontSize',8,'Color','w','YColor','k','XColor','k','XLim',1000*0.05*[-1 1],'YLim',20*[-1 1])
% set(ax([1:4,6]),'XColor','none')
set(ax([2:3,4:6]),'YColor','none')

%% Saccade Velocity %%
FIG = figure (2) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 6 4];
FIG.Name = 'Saccade Velocity';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(N{1,3},1);
for jj = 1:N{1,3}
    ax(jj) = subplot(ceil(N{1,3}/3),3,jj) ; hold on
    ax(jj).Title.String = [num2str(Vel(jj)) ' (' char(176) '/s)'];
    ax(jj).Title.Color = 'k';
    
    cent = max(SACCADE.cIdx{jj},[],'all');
    span = (cent-5):(cent+5);
    
    plot(1000*SACCADE.Head.Time{jj},SACCADE.Head.Velocity{jj}, 'Color', 0.7*CC{jj})
%     plot(1000*SACCADE.HeadStats.Time(jj).Median(span),SACCADE.HeadStats.Velocity(jj).Median(span),'w','LineWidth',3)
    
	PlotPatch(SACCADE.HeadStats.Velocity(jj).Median(span), SACCADE.HeadStats.Velocity(1,1).STD(span), ...
        1000*SACCADE.HeadStats.Time(jj).Median(span), 3, N.fly, CC{jj}, [0.7 0.7 0.7], 0.4, 3);
    
    plot([-100 100],-sign(Vel(jj))*300*[1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',1)
%     plot([-100 100],0*[1 1],'-','Color','k','LineWidth',0.5)
    
    if sign(Vel(jj))==1
        ax(jj).YLim = -1100*[1 -0.2];
    else
        ax(jj).YLim =  1100*[-0.2 1];
    end
    
    if jj~=5
        ax(jj).XTickLabels = '';
        ax(jj).XColor = 'k';
    else
        ax(jj).XLabel.String = 'Time (ms)';
        ax(jj).XColor = 'w';
    end
    
    if (jj-1)/3 ~= round((jj-1)/3)
        ax(jj).YTickLabels = '';
        ax(jj).YColor = 'k';
    else
        ax(jj).YLabel.String = ['(' char(176) '/s)'];
        ax(jj).YColor = 'w';
    end
end
set(ax,'FontSize',8,'Color','w','YColor','k','XColor','k','XLim',1000*0.05*[-1 1])
% set(ax([1:4,6]),'XColor','none')
% set(ax([2:3,5:6]),'YColor','none')

%% Saccade Velocity ALL %%
FIG = figure (20) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 6/3 4/2];
FIG.Name = 'Saccade Velocity';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
clear ax h
ax(1) = subplot(1,1,1) ; hold on
pp = 1;
for jj = 6
    cent = max(SACCADE.cIdx{jj},[],'all');
    span = (cent-5):(cent+5);
    
    h.all{jj} = plot(1000*SACCADE.Head.Time{jj},SACCADE.Head.Velocity{jj}, 'Color', 0.7*CC{jj});
    set(h.all{jj}(:),'Color',0.9*[h.all{jj}(1).Color 0.1])

%     plot(1000*SACCADE.HeadStats.Time(jj).Median(span),SACCADE.HeadStats.Velocity(jj).Median(span),'w','LineWidth',3)
    
	[h.patch(jj),h.grand(jj)] = PlotPatch(SACCADE.HeadStats.Velocity(jj).Median(span), SACCADE.HeadStats.Velocity(1,1).STD(span), ...
        1000*SACCADE.HeadStats.Time(jj).Median(span), 3, N.fly, CC{jj}, [0.7 0.7 0.7], 0.4, 3);
    delete(h.patch(jj))
    
    plot([-100 100],-sign(Vel(jj))*300*[1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',1);
%     plot([-100 100],0*[1 1],'-','Color','k','LineWidth',0.5)
    
    pp = pp + 1;
end
set(ax,'FontSize',8,'Color','w','YColor','k','XColor','k','XLim',1000*0.05*[-1 1],'YLim',[0 820])
ax(1).XLabel.String = 'Time (ms)';
ax.XTick = -50:25:50;
% set(h.all,'Color','c')
% uistack(h.grand(4:6),'top')
% uistack(h.patch(4:6),'top')

%% Inter-Saccade Position %%
FIG = figure (3) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 6 4];
FIG.Name = 'Inter-Saccade Head Position';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gca;
hh = gobjects(N{1,3},1);
hp = gobjects(N{1,3},1);
hold on
SS = repmat({'m','y','c'},1,2);
pLim = 0.70;
for jj = 1:N{1,3}
    % Only plot until a percentage of intervals are still active
    remv = isnan(INTERVAL.Head.Position{jj}(1,:));
    INT = INTERVAL.Head.Position{jj}(:,~remv);
    tLim = sum(isnan(INT),2)./(size(INT,2)-1);
    span = 1:length(tLim(tLim<pLim));
    
    plot(INTERVAL.Head.Time{jj},INTERVAL.Head.Position{jj},'Color', 0.7*CC{jj})
    hh(jj) = plot(INTERVAL.HeadStats.Time(jj).Median(span),INTERVAL.HeadStats.Position(jj).Median(span),...
        'Color', CC{jj}, 'LineWidth', 3);
    
	hp(jj) = plot(INTERVAL.HeadStats.Time(jj).Median(span), Stim(span,jj), '--', 'Color', SS{jj}, 'LineWidth', 2);
    
end
set(ax,'FontSize',8,'Color','w','YColor','k','XColor','k','XLim',[0 1.5])
% ax.YLim = max(abs(ax.YLim))*[-1 1];
ax.YLim = 30*[-1 1];
ax.XLabel.String = 'Time (s)';
ax.YLabel.String = ['Head Angle (' char(176) ')'];

uistack(hh,'top')
uistack(hp,'top')

%% Inter-Saccade Position Error %%
FIG = figure (4) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 3 2];
FIG.Name = 'Inter-Saccade Head Position Error';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gca;
hh = gobjects(N{1,3},1);
hp = gobjects(N{1,3},1);
hold on

SS = repmat({'m','y','c'},1,2);
pLim = 0.7;
for jj = 1:N{1,3}
    % Only plot until a percentage of intervals are still active
    remv = isnan(INTERVAL.Head.Position{jj}(1,:));
    INT = INTERVAL.Head.Position{jj}(:,~remv);
    tLim = sum(isnan(INT),2)./(size(INT,2)-1);
    span = 1:length(tLim(tLim<pLim));
    
    h.trial = plot(INTERVAL.Head.Time{jj},INTERVAL.Head.Position_Error{jj},'Color', 0.5*CC{jj});
    set(h.trial(:),'Color',[h.trial(1).Color 0.2])
    hh(jj) = plot(INTERVAL.HeadStats.Time(jj).Median(span),INTERVAL.HeadStats.Position_Error(jj).Median(span),...
        'Color', CC{jj}, 'LineWidth', 3);
    
%     hp(jj) = plot(INTERVAL.HeadStats.Time(jj).Median(span), Stim(span,jj), '--', 'Color', SS{jj}, 'LineWidth', 2);
    
    plot([0 10],[0 0], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
%     plot([0 ; repmat(INTERVAL.HeadStats.Time(jj).Median(span(end)),2,1)],...
%         [0 , INTERVAL.HeadStats.Position_Error(jj).Median(span(end)) , 0],...
%         '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
    
	plot([repmat(INTERVAL.HeadStats.Time(jj).Median(span(end)),2,1)],...
        [INTERVAL.HeadStats.Position_Error(jj).Median(span(end)) , 0],...
        '--', 'Color', CC{jj}, 'LineWidth', 2);
    
end
set(ax,'FontSize',8,'Color','w','YColor','k','XColor','k','XLim',[0 1.5])
% ax.YLim = max(abs(ax.YLim))*[-1 1];
ax.YLim = 60*[-1 1];
ax.XLabel.String = 'Time (s)';
ax.YLabel.String = ['Head Position Error (' char(176) ')'];

uistack(hh,'top')
uistack(hp,'top')

%% Inter-Saccade Velocity %%
FIG = figure (5) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 6 4];
FIG.Name = 'Inter-Saccade Head Velocity';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gca;
hh = gobjects(N{1,3},1);
hp = gobjects(N{1,3},1);
hold on
pLim = 0.90;
for jj = 1:N{1,3}
	ax(jj) = subplot(ceil(N{1,3}/3),3,jj) ; hold on
    ax(jj).Title.String = [num2str(Vel(jj)) ' (' char(176) '/s)'];
    ax(jj).Title.Color = 'k';
    
    % Only plot until a percentage of intervals are still active
    remv = isnan(INTERVAL.Head.Position{jj}(1,:));
    INT = INTERVAL.Head.Position{jj}(:,~remv);
    tLim = sum(isnan(INT),2)./(size(INT,2)-1);
    span = 1:length(tLim(tLim<pLim));
    
%     plot(INTERVAL.Head.Time{jj},INTERVAL.Head.Velocity{jj},'Color', 0.5*CC{jj})
    hh(jj) = plot(INTERVAL.HeadStats.Time(jj).Median(span),INTERVAL.HeadStats.Velocity(jj).Median(span),...
        'Color', CC{jj}, 'LineWidth', 2);
    
    hp(jj) = plot([0 10], repmat(Vel(jj),1,2), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
    
    ax(jj).XLabel.String = 'Time (s)';
    ax(jj).YLabel.String = ['Head Velocity (' char(176) '/s)'];
    
    ax(jj).YLim = 150*[-1 1];
    
end
set(ax,'FontSize',8,'Color','w','YColor','k','XColor','k','XLim',[0 0.5])

%% Inter-Saccade Velocity Error %%
FIG = figure (6) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 6 4];
FIG.Name = 'Inter-Saccade Head Velocity Error';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gca;
hh = gobjects(N{1,3},1);
hp = gobjects(N{1,3},1);
hold on
pLim = 0.90;
for jj = 1:N{1,3}
	ax(jj) = subplot(ceil(N{1,3}/3),3,jj) ; hold on
    ax(jj).Title.String = [num2str(Vel(jj)) ' (' char(176) '/s)'];
    ax(jj).Title.Color = 'k';
    
    % Only plot until a percentage of intervals are still active
    remv = isnan(INTERVAL.Head.Position{jj}(1,:));
    INT = INTERVAL.Head.Position{jj}(:,~remv);
    tLim = sum(isnan(INT),2)./(size(INT,2)-1);
    span = 1:length(tLim(tLim<pLim));
    
%     plot(INTERVAL.Head.Time{jj},INTERVAL.Head.Velocity_Error{jj},'Color', 0.5*CC{jj})
    hh(jj) = plot(INTERVAL.HeadStats.Time(jj).Median(span),INTERVAL.HeadStats.Velocity_Error(jj).Median(span),...
        'Color', CC{jj}, 'LineWidth', 2);
    
    hp(jj) = plot([0 10], repmat(Vel(jj),1,2), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
    
    ax(jj).XLabel.String = 'Time (s)';
    ax(jj).YLabel.String = ['Head Velocity Error (' char(176) '/s)'];
    
	if sign(Vel(jj))==1
        ax(jj).YLim = -100*[0.2 -1];
    else
        ax(jj).YLim =  100*[-1 0.2];
    end
    
end
set(ax,'FontSize',8,'Color','w','YColor','k','XColor','k','XLim',[0 1.5])

%% Saccade Trigger Polar Plot %%
FIG = figure (7) ; clf
FIG.Units = 'inches';
FIG.Position = 0.75*[2 2 2 2];
FIG.Name = 'Saccade Trigger';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';

polarhistogram([deg2rad(SACD.Head.StartPos(SACD.Head.Dir==-1));-deg2rad(SACD.Head.StartPos(SACD.Head.Dir==1))],100,'FaceColor','g','FaceAlpha',.9); hold on
polarhistogram([deg2rad(SACD.Head.EndPos(SACD.Head.Dir==-1));-deg2rad(SACD.Head.EndPos(SACD.Head.Dir==1))],100,'FaceColor','r','FaceAlpha',.9) ; hold on
% polarhistogram(-deg2rad(SACD.Head.StartPos(SACD.Head.Dir==1)),100,'FaceColor','g','FaceAlpha',.9); hold on
% polarhistogram(-deg2rad(SACD.Head.EndPos(SACD.Head.Dir==1)),100,'FaceColor','r','FaceAlpha',.9) ; hold on
% polarhistogram(deg2rad(SACD.Head.StartPos),100,'FaceColor','g','FaceAlpha',.9); hold on
% polarhistogram(deg2rad(SACD.Head.EndPos),100,'FaceColor','r','FaceAlpha',.9) ; hold on
ax = gca;
ax.FontSize = 8;
grid off
axis tight
ax.Color = 'w';
thetalim([-20 20])
thetaticks([-20 -10 0 10 20])
thetaticklabels({'20','10','0','-10','-20'})
% d = ax.ThetaDir;
ax.ThetaZeroLocation = 'top';
ax.ThetaAxis.Label.String = ['Head Position (' char(176) ')'];
ax.ThetaAxis.Label.FontSize = 8;
% leg = legend('Start','End','Location','North');
% leg.Location = 'northwest';
% leg.Box = 'off';

%% Saccade Statistics %%
G = [];
pp = 1:3;
for jj = 1:3
    [rr,~] = find( SACD.Head.speed == U{1,3}{1}(jj) );
    G(rr,1) = jj;
end

SS = [5,6,14,16,18,23,25,24,26,22];
YY = {  ['Duration (ms)'],...
        ['Amplitude (' char(176) ')'],...
        ['Trigger Position (' char(176) ')'],...
        ['End Position (' char(176) ')'],...
        ['Peak Velocity (' char(176) '/s)'],...
        ['Error (' char(176) ')'],...
        ['Integrated Error (' char(176) ' \ast s)'],...
      	['Velocity Err (' char(176) '/s)'],...
        ['Integrated Velocity Err (' char(176) ')'],...
        ['Rate (#/s)'] };
    
CC = {[0 0 0.5],[0 0 0],[0 0.7 0 ],[0.7 0 0],[0.5 0.5 0.5],[0.1 0.5 0.5],[0.5 0.1 0.5],[0.9 0.7 0.7],[0.3 0.3 0.9],...
    [0.6 0.1 0.4],[0 0.8 0.8],[0.9 0.7 0.3]};
FF = gobjects(length(SS),1);
AX = gobjects(length(SS),1);
for ww = 1:length(SS)
    FF(ww) = figure (10+ww) ; clf ; hold on
    FF(ww).Units = 'inches';
    FF(ww).Position = [0+(0.6*ww) , 6-(0.3*ww) , 2 , 2];
    FF(ww).Color = 'w';
    AX(ww) = gca;
    axis tight
        data = SACD.Head{:,SS(ww)};
        if any(SS(ww)==[5,6,14,16,18,23,25,24,26,22])
            data = abs(data);
        end
        bx = boxplot(data,G,'Labels',{Vel(pp)},'Width',0.5,'Symbol','','Whisker',2);
        xlabel(['Stimulus Velocity (' char(176) '/s)'])
        ylabel(YY{ww})
        
        h = get(bx(5,:),{'XData','YData'});
        for k = 1:size(h,1)
           patch(h{k,1},h{k,2},CC{ww});
        end
        
        set(findobj(FF(ww),'tag','Median'), 'Color', 'w');
        set(findobj(FF(ww),'tag','Box'), 'Color', 'k');
        set(findobj(FF(ww),'tag','Upper Whisker'), 'Color', 'k');
        AX(ww).Children = AX(ww).Children([end 1:end-1]);
        
        uLim = findobj(FF(ww),'tag','Upper Whisker');
        lLim = findobj(FF(ww),'tag','Lower Whisker');
        uLim = uLim(3).YData(2);
        lLim = lLim(1).YData(1);
        AX(ww).YLim = [0.5*lLim 1.5*uLim];
        
end
set(AX,'FontSize',8)

%% Saccade Removal %%
FIG = figure (20) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 6 4];
FIG.Name = 'Removed Saccades';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
hold on
clear h
ax = gobjects(N{1,3},1);
rmvIdx = 4;
for jj = 1:N.vel
    ax(jj) = subplot(2,3,jj) ; hold on
    ax(jj).Title.String = ['Velocity (' char(176) '/s)'];
    ax(jj).XLabel.String = 'Time(s)';
    ax(jj).YLabel.String = ['Head Angle (' char(176) ')'];
    for kk = 1:N.fly
        for ii = 1:size(TRIAL{kk,jj},1)
            h.trial = plot(TRIAL{kk,jj}{ii,rmvIdx}.Time, TRIAL{kk,jj}{ii,rmvIdx}.X(:,1), 'Color', 0.5*CC{jj});
            h.trial.Color(4) = 0.2;
        end
        plot(FLY{jj}{kk,rmvIdx}.Mean{5}, FLY{jj}{kk,rmvIdx}.Mean{6}(:,1), 'Color', 0.7*CC{jj},'LineWidth',1)
    end
    plot(GRAND{jj,rmvIdx}.Mean{1}{5}, GRAND{jj,rmvIdx}.Mean{1}{6}(:,1), 'Color', CC{jj},'LineWidth',3)
    plot(GRAND{jj,rmvIdx}.Mean{1}{5},Stim(:,jj), '--', 'Color', 'k', 'LineWidth',2)
end
set(ax,'XLim',[0 10])
set(ax(1:3),'YLim',200*[-0.1 1])
set(ax(4:6),'YLim',200*[-1 0.1])
set(ax([1:4,6]),'XColor','none')
set(ax([2:3,5:6]),'YColor','none')
set(ax,'FontSize',8)

%% Saccde Frequency Domain %%
FIG = figure (21) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 6 4];
FIG.Name = 'Full Trial Frequency Domain';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(N{1,3},1);
for jj = 1:N{1,3}
    ax(jj) = subplot(2,3,jj); hold on
    for kk = 1:N.fly
        for ii = 1:size(TRIAL{kk,jj},1)
            plot(TRIAL{kk,jj}{ii,2}.Fv, TRIAL{kk,jj}{ii,2}.Mag(:,1), 'Color', 0.5*CC{jj})
        end
        plot(FLY{jj}{kk,2}.Mean{7}, FLY{jj}{kk,2}.Mean{8}(:,1), 'Color', CC{jj},'LineWidth',1)
    end
    plot(GRAND{jj,2}.Mean{1}{7}, GRAND{jj,2}.Mean{1}{8}(:,1), 'Color', 'k','LineWidth',2)
end
set(ax,'XLim',[0 40])
set(ax,'YLim',[0 1.2])



%% Hypothesis Testing
test1 = abs(SACD.Head.PeakVel(SACD.Head.speed==8));
test2 = abs(SACD.Head.PeakVel(SACD.Head.speed==16));

clc
[h,p] = ttest2(test1,test2)


end