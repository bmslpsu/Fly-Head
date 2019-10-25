function [FIG] = MakeFig_Static_Head_Saccade_Stats()
%% MakeFig_Static_Head_Saccade_Stats:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

load(fullfile(root,FILE{1}),'SACD','U','N','T');

clearvars -except SACD U N T Wave

Wave = U.SpatFreq{1};

%% Saccade Statistics %%
pp = 1:N.SpatFreq;
n_sacd = size(SACD.Head,1);
G = nan(n_sacd,1);
for jj = pp
    waves = SACD.Head.WaveIdx;    
    [rr,~] = find( waves == pp(jj) );
    G(rr,1) = jj;
end

SS = [6,7,19,21,15,17,23];
n_plot = length(SS);
clms = 3;
rows = ceil(n_plot/clms);

YY = {  ['Duration (ms)'],...
        ['Amplitude (' char(176) ')'],...
        ['Peak Velocity (' char(176) '/s)'],...
        ['Threshold (' char(176) '/s)'],...
        ['Trigger Position (' char(176) ')'],...
        ['End Position (' char(176) ')'],...
        ['Rate (#/s)']};
    
CC = {[0 0 0.5],[0 0 0],[0 0.7 0 ],[0.7 0 0],[0.5 0.5 0.5],[0.1 0.5 0.5],...
      [0.5 0.1 0.5],[0.9 0.7 0.7],[0.3 0.3 0.9],[0.6 0.1 0.4],[0 0.8 0.8],[0.9 0.7 0.3]};

FIG = figure (1) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 2.5*clms 2*rows];
movegui(FIG,'center')
ax = gobjects(length(SS),1);
for ww = 1:n_plot-1
    ax(ww) = subplot(rows,clms,ww); axis tight
    data = SACD.Head{:,SS(ww)};
    if any(SS(ww)==[6,7,19,21,15,17,23])
        data = abs(data);
    end
    if ww==1
        data = 1000*data;
    end
    
    bx = boxplot(data,G,'Labels',{Wave(pp)},'Width',0.5,'Symbol','','Whisker',2);
    xlabel(['Stimulus Velocity (' char(176) '/s)'])
    ylabel(YY{ww})
    box off
%     ax(ww).XColor = 'none';
%     ax(ww).XAxis.Label.Color = 'k';
%     ax(ww).XAxis.Label.Visible = 'on';

    h = get(bx(5,:),{'XData','YData'});
    for k = 1:size(h,1)
       patch(h{k,1},h{k,2},CC{ww});
    end

    set(findobj(ax(ww),'tag','Median'), 'Color', 'w','LineWidth',1.5);
    set(findobj(ax(ww),'tag','Box'), 'Color', 'none');
    set(findobj(ax(ww),'tag','Upper Whisker'), 'Color', 'k');
    ax(ww).Children = ax(ww).Children([end 1:end-1]);

    uLim = findobj(ax(ww),'tag','Upper Whisker');
    lLim = findobj(ax(ww),'tag','Lower Whisker');
    uLim = uLim(3).YData(2);
    lLim = lLim(1).YData(1);
    ax(ww).YLim = [0.5*lLim 1.5*uLim];
    ax(ww).YLim(1) = 0;
end
ax(1).YLim(2) = 0.08e3;
ax(2).YLim(2) = 30;
ax(3).YLim(2) = 1000;
ax(4).YLim(2) = 600;
ax(5).YLim(2) = 20;
ax(6).YLim(2) = 20;

ax(n_plot) = subplot(rows,clms,n_plot); axis tight
data = SACD.Head{:,SS(n_plot)};  
bx = boxplot(data,G,'Labels',{Wave(pp)},'Width',0.5,'Symbol','','Whisker',2);
xlabel(['Stimulus Velocity (' char(176) '/s)'])
ylabel(YY{n_plot})
box off

h = get(bx(5,:),{'XData','YData'});
for k = 1:size(h,1)
   patch(h{k,1},h{k,2},CC{n_plot});
end

set(findobj(ax(n_plot),'tag','Median'), 'Color', 'w','LineWidth',1.5);
set(findobj(ax(n_plot),'tag','Box'), 'Color', 'none');
set(findobj(ax(n_plot),'tag','Upper Whisker'), 'Color', 'k');
ax(n_plot).Children = ax(n_plot).Children([end 1:end-1]);

uLim = findobj(ax(n_plot),'tag','Upper Whisker');
lLim = findobj(ax(n_plot),'tag','Lower Whisker');
uLim = uLim(3).YData(2);
lLim = lLim(1).YData(1);
ax(n_plot).YLim = [0.5*lLim 1.5*uLim];
ax(n_plot).YLim(1) = 0;
ax(n_plot).YLim(2) = 2;

set(ax(1:6),'FontSize',8)

end