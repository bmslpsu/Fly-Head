function [FIG] = MakeFig_SOS_Walking_head_Time_all()
%% MakeFig_SOS_HeadFree_pat_head_Time:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'E:\Walking_Experiments';

[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

Walking = load(fullfile(root,FILE{1}),'TRIAL','GRAND', 'FLY','U','N', 'T');
%%
figNum = 1;
catIdx = 2;
xIdx = 1;

filename = 'SOS_Walking_all_head_Time';

FIG = figure (figNum); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 8 2];
movegui(FIG,'center')
FIG.Name = filename;

ax.L = gca; hold on
ax.L.YColor = 'k';
ax.L.FontSize = 12;
ax.L.YLabel.String = ['Head (' char(176) ')'];
ax.L.YLabel.FontSize = 12;
ax.L.YLim = 20*[-1 1];
ax.L.YTick = 20*[-1 0 1];
ax.L.XLabel.String = 'Time (s)';
ax.L.XLabel.Color = 'k';
ax.L.XLabel.FontSize = 12;
ax.L.XLim = [0 20];

% for kk = 1:Walking.N{1,1}
%     for jj = 1:Walking.T{kk,2}
%     plot(Walking.TRIAL{kk}{jj,catIdx}.Time,Walking.TRIAL{kk}{jj,catIdx}.X(:,xIdx),'color',[.5 .5 .5 .5], 'LineWidth',.25)
%     end
% end
% 
% for kk = 1:Walking.N{1,1}
%      plot(Walking.FLY{kk,catIdx}.Mean{5},Walking.FLY{kk,catIdx}.Mean{6}(:,1),'LineWidth',.5)
% end
 plot(Walking.TRIAL{1}{1,1}.Time,Walking.TRIAL{1}{1,1}.X(:,xIdx),'color', [0 0.4470 0.7410],'LineWidth',1)
 plot(Walking.GRAND{1,catIdx}.Mean{2}{5},Walking.GRAND{1,catIdx}.Mean{2}{6}(:,xIdx),'k','LineWidth',2)
 legend('Pattern', 'Fly Grand Mean')
 

end