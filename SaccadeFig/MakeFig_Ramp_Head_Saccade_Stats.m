function [FIG] = MakeFig_Ramp_Head_Saccade_Stats()
%% MakeFig_Ramp_Head_Saccade_Stats:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG    :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

[FILE.Anti,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE.Anti = cellstr(FILE.Anti)';

% [FILE.Static,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
%     'Select head angle trials', root, 'MultiSelect','off');
% FILE.Static = cellstr(FILE.Static)';

Anti    = load(fullfile(root,FILE.Anti{1}),'SACD','U','N','T');
% Static  = load(fullfile(root,FILE.Anti{1}),'SACD','U','N','T');

Vel = Anti.U{1,3}{1};

%% Saccade Statistics %%
match = -1;
mIdx = Anti.SACD.Head.Match==match;
% mIdx = (Anti.SACD.Head.Match==1) | (Anti.SACD.Head.Match==-1); % for anti & co -directional
pp = 1:(length(Anti.U{1,3}{1})/2); % for each speed
% pp = 1:N.vel; % for each velocity
w_scale = 2*length(pp)/10;
n_sacd = sum(mIdx);
G = nan(n_sacd,1);
speeds = Anti.SACD.Head.velIdx(mIdx);
for jj = pp
    [rr,~] = find( speeds == pp(jj) );
    G(rr,1) = jj;
end

SS = [6,7,19,24,26,15,17,25,27,21,23];
n_plot = length(SS);
clms = 3;
rows = ceil(n_plot/clms);

YY = {  ['Duration (ms)'],...
        ['Amplitude (' char(176) ')'],...
        ['Peak Velocity (' char(176) '/s)'],...
        ['Error (' char(176) ')'],...
        ['Integrated Error (' char(176) ' \ast s)'],...
        ['Trigger Position (' char(176) ')'],...
        ['End Position (' char(176) ')'],...
      	['Velocity Err (' char(176) '/s)'],...
        ['Integrated Velocity Err (' char(176) ')'],...
        ['Threshold (' char(176) '/s)'],...
        ['Rate (#/s)']};
    
% CC = {[0 0 0.5],[0 0 0],[0 0.7 0 ],[0.7 0 0],[0.5 0.5 0.5],[0.1 0.5 0.5],...
%       [0.5 0.1 0.5],[0.9 0.7 0.7],[0.3 0.3 0.9],[0.6 0.1 0.4],[0 0.8 0.8],[0.9 0.7 0.3]};
CC = repmat(hsv(Anti.N.vel/2),2,1);

FIG = figure (1) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 2*clms*w_scale 2*rows];
movegui(FIG,'center')
ax = gobjects(length(SS),1);
for ww = 1:n_plot-2
    ax(ww) = subplot(rows,clms,ww); axis tight
    data = Anti.SACD.Head{mIdx,SS(ww)};
    if any(SS(ww)==[6,7,19,24,26,25,27,21,23])
        data = abs(data);
    end
    if ww==1
        data = 1000*data;
    end
    
    bx = boxplot(data,G,'Labels',{Vel(pp)},'Width',0.5,'Symbol','','Whisker',2);
    xlabel(['Stimulus Velocity (' char(176) '/s)'])
    ylabel(YY{ww})
    box off

    h = get(bx(5,:),{'XData','YData'});
    for kk = 1:size(h,1)
       patch(h{kk,1},h{kk,2},CC(kk,:));
    end

    set(findobj(ax(ww),'tag','Median'), 'Color', 'w','LineWidth',1.5);
    set(findobj(ax(ww),'tag','Box'), 'Color', 'none');
    set(findobj(ax(ww),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
    set(findobj(ax(ww),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
    ax(ww).Children = ax(ww).Children([end 1:end-1]);
    ax(ww).YLim(1) = 0;
end


fIdx = find(~isnan(Anti.SACD.Head.Rate(mIdx)));

for ww = (n_plot-1):n_plot
    data = Anti.SACD.Head{mIdx,SS(ww)};
    data = data(fIdx);
    ax(ww) = subplot(rows,clms,ww); axis tight
    bx = boxplot(data,G(fIdx),'Labels',{Vel(pp)},'Width',0.5,'Symbol','','Whisker',2);
    xlabel(['Stimulus Velocity (' char(176) '/s)'])
    ylabel(YY{ww})
    box off
    
    h = get(bx(5,:),{'XData','YData'});
    for kk = 1:size(h,1)
       patch(h{kk,1},h{kk,2},CC(kk,:));
    end
    
    set(findobj(ax(ww),'tag','Median'), 'Color', 'w','LineWidth',1.5);
    set(findobj(ax(ww),'tag','Box'), 'Color', 'none');
    set(findobj(ax(ww),'tag','Upper Whisker'), 'Color', 'k');
    ax(ww).Children = ax(ww).Children([end 1:end-1]);

end

set(ax,'FontSize',8)

ax(1).YLim(2)   = 0.08e3;
ax(2).YLim(2)   = 30;
ax(3).YLim(2)   = 1000;
ax(4).YLim(2)   = 160;
ax(5).YLim(2)   = 70;
ax(6).YLim      = 20*[-1 1];
ax(7).YLim      = ax(6).YLim;
ax(8).YLim(2)   = 350;
ax(9).YLim(2)   = 500;
ax(10).YLim(2)  = 700;
ax(11).YLim(2)  = 2.2;

end