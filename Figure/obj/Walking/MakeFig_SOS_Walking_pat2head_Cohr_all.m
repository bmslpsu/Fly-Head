function [FIG] = MakeFig_SOS_Walking_pat2head_Cohr_all()
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
figNum = 2;
catIdx = 4;
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
ax.L.FontSize = 12 ;
ax.L.YLabel.String = ['Magnitude'];
ax.L.YLabel.FontSize = 12;
ax.L.XLabel.String = 'Frequency';
ax.L.YLim = [0 1];
ax.L.XLabel.Color = 'k';
ax.L.XLabel.FontSize = 12;
ax.L.XLim = [0 12];
ax.L.XTick = [1, 3.1, 5.3, 7.4, 9.6];

% for kk = 1:Walking.N{1,1}
%     for jj = 1:Walking.T{kk,2}
%     plot(Walking.TRIAL{kk}{jj,catIdx}.CoherenceFV,Walking.TRIAL{kk}{jj,catIdx}.Coherence,'color',[.5 .5 .5 .5], 'LineWidth',.5)
%     end
% end
% 
% for kk = 1:Walking.N{1,1}
%      plot(Walking.FLY{kk,catIdx}.Mean{8},Walking.FLY{kk,catIdx}.Mean{7}(:,xIdx),'LineWidth',1)
% end

k.plot = PlotPatch(Walking.GRAND{1,catIdx}.Mean{2}{9}(:,xIdx), Walking.GRAND{1, catIdx}.STD{2}{9}(:, xIdx), ...
    Walking.GRAND{1,catIdx}.Mean{2}{4}, 2, Walking.N{1,1},[0 0.4470 0.7410],[0 0.4470 0.7410], 0.5, 3);
h.plot = PlotPatch(Walking.GRAND{1,catIdx}.Mean{2}{7}(:,xIdx),Walking.GRAND{1,catIdx}.STD{2}{7}(:,xIdx), ...
    Walking.GRAND{1,catIdx}.Mean{2}{8}(:,xIdx),2,Walking.N{1,1},'k',[0.4 0.4 0.6],0.5,3);


end