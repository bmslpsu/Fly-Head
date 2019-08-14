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

load(fullfile(root,FILE{1}),'SACCADE','INTERVAL','SACD','Stim','U','N','I','TRIAL');

CC = repmat({[1 0 0],[0 1 0],[0 0 1]},1,2);

Vel = 3.75*U{1,3}{1};

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
set(ax([1:4,6]),'XColor','none')
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
        ax(jj).YLim = -1100*[1 -1];
    else
        ax(jj).YLim =  1100*[-1 1];
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
set(ax([1:4,6]),'XColor','none')
set(ax([2:3,5:6]),'YColor','none')

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
    
	hp(jj) = plot(INTERVAL.HeadStats.Time(jj).Median(span), Stim(span,jj), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
    
end
set(ax,'FontSize',8,'Color','w','YColor','k','XColor','k','XLim',[0 2])
% ax.YLim = max(abs(ax.YLim))*[-1 1];
ax.YLim = 30*[-1 1];
ax.XLabel.String = 'Time (s)';
ax.YLabel.String = ['Head Angle (' char(176) ')'];

uistack(hh,'top')
uistack(hp,'top')

%% Inter-Saccade Position Error %%
FIG = figure (4) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 6 4];
FIG.Name = 'Inter-Saccade Head Position Error';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gca;
hh = gobjects(N{1,3},1);
hp = gobjects(N{1,3},1);
hold on
pLim = 0.7;
for jj = 1:N{1,3}
    % Only plot until a percentage of intervals are still active
    remv = isnan(INTERVAL.Head.Position{jj}(1,:));
    INT = INTERVAL.Head.Position{jj}(:,~remv);
    tLim = sum(isnan(INT),2)./(size(INT,2)-1);
    span = 1:length(tLim(tLim<pLim));
    
    plot(INTERVAL.Head.Time{jj},INTERVAL.Head.Position_Error{jj},'Color', 0.5*CC{jj})
    hh(jj) = plot(INTERVAL.HeadStats.Time(jj).Median(span),INTERVAL.HeadStats.Position_Error(jj).Median(span),...
        'Color', CC{jj}, 'LineWidth', 3);
    
    hp(jj) = plot(INTERVAL.HeadStats.Time(jj).Median(span), Stim(span,jj), '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
    
end
set(ax,'FontSize',8,'Color','w','YColor','k','XColor','k','XLim',[0 2])
ax.YLim = max(abs(ax.YLim))*[-1 1];
% ax.YLim = 28*[-1 1];
ax.XLabel.String = 'Time (s)';
ax.YLabel.String = ['Head Angle (' char(176) ')'];

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
    ax(jj).YLabel.String = ['Head Velocity (' char(176) '/s)'];
    
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
FIG.Position = [2 2 4 4];
FIG.Name = 'Saccade Trigger';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';

polarhistogram(deg2rad(SACD.Head.StartPos(SACD.Head.Dir==1)),100,'FaceColor','g','FaceAlpha',.9); hold on
polarhistogram(deg2rad(SACD.Head.EndPos(SACD.Head.Dir==1)),100,'FaceColor','r','FaceAlpha',.9) ; hold on
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
leg = legend('Start','End','Location','North');
leg.Location = 'northwest';

%% Saccade Statistics %%
G = [];
pp = 1:3;
for jj = 1:3
    [rr,~] = find( SACD.Head.speed == U{1,3}{1}(jj) );
    G(rr,1) = jj;
end

SS = [5,6,14,16,21,25,27,26,28];
YY = {  ['Duration (ms)'],...
        ['Amplitude (' char(176) ')'],...
        ['Trigger Position (' char(176) ')'],...
        ['End Position (' char(176) ')'],...
        ['Peak Velocity (' char(176) '/s)'],...
        ['Error (' char(176) ')'],...
        ['Integrated Error (' char(176) ' \ast s)'],...
      	['Velocity Err (' char(176) '/s)'],...
        ['Integrated Velocity Err (' char(176) ')'],   };
    
CC = {[0 0 0.5],[0 0 0],[0 0.7 0 ],[0.7 0 0],[0.5 0.5 0.5],[0.1 0.5 0.5],[0.5 0.1 0.5],[0.9 0.7 0.7],[0.3 0.3 0.9],...
    [0.6 0.1 0.8],[0 0.8 0.8]};
FF = gobjects(length(SS),1);
AX = gobjects(length(SS),1);
for ww = 1:length(SS)
    FF(ww) = figure (10+ww) ; clf ; hold on
    FF(ww).Units = 'inches';
    FF(ww).Position = [0+(0.6*ww) , 6-(0.3*ww) , 3 , 3];
    FF(ww).Color = 'w';
    AX(ww) = gca;
    axis tight
        data = SACD.Head{:,SS(ww)};
        if any(SS(ww)==[6,14,16,25,26,27,28])
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

%% Hypothesis Testing
test1 = abs(SACD.Head.PeakVel(SACD.Head.speed==8));
test2 = abs(SACD.Head.PeakVel(SACD.Head.speed==16));

clc
[h,p] = ttest2(test1,test2)


end