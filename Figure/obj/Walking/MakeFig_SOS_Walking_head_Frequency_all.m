function [FIG] = MakeFig_SOS_Walking_head_Frequency_all()
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
ax.L.YLabel.String = ['Magnitude'];
ax.L.YLabel.FontSize = 12;
% ax.L.YLim =[0 10];
% ax.L.YTick = 15*[-1 0 1];
ax.L.XLabel.String = 'Frequency';
ax.L.XLabel.Color = 'k';
ax.L.XLabel.FontSize = 12;
ax.L.XLim = [0 20];
ax.L.XTick = [1, 3.1, 5.3, 7.4, 9.6, 20];

for kk = 1:Walking.N{1,1}
    for jj = 1:Walking.T{kk,2}
    plot(Walking.TRIAL{kk}{jj,catIdx}.Fv,Walking.TRIAL{kk}{jj,catIdx}.Mag(:,xIdx),'color',[.5 .5 .5 .5], 'LineWidth',.5)
    end
end

for kk = 1:Walking.N{1,1}
     plot(Walking.FLY{kk,catIdx}.Mean{7},Walking.FLY{kk,catIdx}.Mean{8}(:,xIdx),'LineWidth',1)
end
plot(Walking.TRIAL{1}{1,1}.IOFreq, Walking.TRIAL{1}{1,1}.IOMag(:,xIdx),'color', [0.6350 0.0780 0.1840])
plot(Walking.TRIAL{1}{1,1}.Fv,Walking.TRIAL{1}{1,1}.Mag(:,xIdx),'color', [0.6350 0.0780 0.1840],'LineWidth',2)
plot(Walking.GRAND{1,catIdx}.Mean{2}{7},Walking.GRAND{1,catIdx}.Mean{2}{8}(:,xIdx),'k','LineWidth',2)
plot(Walking.GRAND{1,catIdx}.Mean{2}{10},Walking.GRAND{1,catIdx}.Mean{2}{11}(:,xIdx),'k')
end