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

load(fullfile(root,FILE{1}),'SACCADE','INTERVAL','SACD','Stim','U','N');

filename = 'Ramp_Head_Saccade';

CC = repmat({'r','g','b'},1,2);

%% Saccade Position %%
FIG = figure (1) ; clf
FIG.Position = [100 100 1200 800]*0.75;
FIG.Name = filename;
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'k';
ax = axes;
for jj = 1:N{1,3}
    Vel = 3.75*U{1,3}{1}(jj);
    ax(jj) = subplot(ceil(N{1,3}/3),3,jj) ; hold on
    ax(jj).Title.String = [num2str(Vel) ' (' char(176) '/s)'];
    ax(jj).Title.Color = 'w';
    
    cent = max(SACCADE.cIdx{jj},[],'all');
    span = (cent-5):(cent+5);
    
    plot(1000*SACCADE.Head.Time{jj},SACCADE.Head.Position{jj},CC{jj})
    plot(1000*SACCADE.HeadStats.Time(jj).Median(span),SACCADE.HeadStats.Position(jj).Median(span),'w','LineWidth',3)
    
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
set(ax,'FontSize',12,'Color','k','YColor','w','XColor','w','XLim',1000*0.05*[-1 1],'YLim',20*[-1 1])

%% Saccade Velocity %%
FIG = figure (2) ; clf
FIG.Position = [100 100 1200 800]*0.75;
FIG.Name = filename;
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'k';
ax = axes;
for jj = 1:N{1,3}
    ax(jj) = subplot(ceil(N{1,3}/3),3,jj) ; hold on
    Vel = 3.75*U{1,3}{1}(jj);
    ax(jj).Title.String = [num2str(Vel) ' (' char(176) '/s)'];
    ax(jj).Title.Color = 'w';
    
    cent = max(SACCADE.cIdx{jj},[],'all');
    span = (cent-5):(cent+5);
    
    plot(1000*SACCADE.Head.Time{jj},SACCADE.Head.Velocity{jj},CC{jj})
    plot(1000*SACCADE.HeadStats.Time(jj).Median(span),SACCADE.HeadStats.Velocity(jj).Median(span),'w','LineWidth',3)
%     plot([-100 100],-sign(Vel)*500*[1 1],'--g','LineWidth',2)
    
    if sign(Vel)==1
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
set(ax,'FontSize',12,'Color','k','YColor','w','XColor','w','XLim',1000*0.05*[-1 1])

%% Inter-Saccade Position %%
FIG = figure (3) ; clf
FIG.Position = [100 100 1200 800]*0.75;
FIG.Name = filename;
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'k';
ax = axes;
hold on
for jj = 1:N{1,3}
    Vel = 3.75*U{1,3}{1}(jj);
%     ax(jj) = subplot(ceil(N{1,3}/3),3,jj) ; hold on
%     ax(jj).Title.String = [num2str(Vel) ' (' char(176) '/s)'];
%     ax(jj).Title.Color = 'w';
    
    span = 1:500;
    
%     plot(INTERVAL.Head.Time{jj},INTERVAL.Head.Velocity{jj},'Color',CC{jj})
    hh(jj) = plot(INTERVAL.HeadStats.Time(jj).Median(span),INTERVAL.HeadStats.Velocity(jj).Median(span),CC{jj},'LineWidth',1)
    
%     if jj<=N{1,3}/2
%         ax(jj).XTickLabels = '';
%     else
%         ax(jj).XLabel.String = 'Time (ms)';
%     end
%     
%     if (jj-1)/3 ~= round((jj-1)/3)
%         ax(jj).YTickLabels = '';
%     else
%         ax(jj).YLabel.String = ['(' char(176) ')'];
%     end
end
set(ax,'FontSize',12,'Color','k','YColor','w','XColor','w','XLim',[0 5])

% ,'YLim',30*[-1 1]

uistack(hh,'top')






%% Saccade Statistics %%
G = [];
pp = 1:3;
for jj = 1:3
    [rr,~] = find( SACD.Head.speed == U{1,3}{1}(jj) );
    G(rr,1) = jj;
end

SS = [5,6,14,16,21,25,27,26,28,5];
YY = {  ['Duration (ms)'],...
        ['Amplitude (' char(176) ')'],...
        ['Trigger Position (' char(176) ')'],...
        ['End Position (' char(176) ')'],...
        ['Peak Velocity (' char(176) '/s)'],...
        ['Error (' char(176) ')'],...
        ['Integrated Error (' char(176) ' \ast s)'],...
      	['Velocity Err (' char(176) '/s)'],...
        ['Integrated Velocity Err (' char(176) ')'],...
        ['Inter-Saccade Duration (ms)']     };
    
CC = {[0 0 0.5],[0 0 0],[0 0.7 0 ],[0.7 0 0],[0.5 0.5 0.5],[0.1 0.5 0.5],[0.5 0.1 0.5],[0.9 0.7 0.7],[0.3 0.3 0.9],...
    [0.6 0.1 0.8],[0 0.8 0.8]};
FF = cell(length(SS),1);
AX = cell(length(SS),1);
for ww = 1:length(SS)
    FF{ww} = figure (2+ww) ; clf ; hold on
    FF{ww}.Position = [-100+(100*ww) , 500-(20*ww) , 500 , 250];
    FF{ww}.Color = 'w';
    AX{ww} = gca;
    axis tight
        data = SACD.Head{:,SS(ww)};
        if SS(ww)==6
            data = abs(data);
        end
        bx = boxplot(data,G,'Labels',{num2str(U{1,3}{1}(pp))},'Width',0.5,'Symbol','','Whisker',2);
        ax = gca;
        xlabel(['Stimulus Velocity (' char(176) '/s)'])
        ylabel(YY{ww})
        set(gca,'FontSize',12)
        h = get(bx(5,:),{'XData','YData'});
        for k=1:size(h,1)
           patch(h{k,1},h{k,2},CC{ww});
        end
        set(findobj(gcf,'tag','Median'), 'Color', 'w');
        set(findobj(gcf,'tag','Box'), 'Color', 'k');
        set(findobj(gcf,'tag','Upper Whisker'), 'Color', 'k');
        ax.Children = ax.Children([end 1:end-1]);
end
end