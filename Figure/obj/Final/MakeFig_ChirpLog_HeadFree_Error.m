function [FIG] = MakeFig_ChirpLog_HeadFree_Error()
%% MakeFig_ChirpLog_HeadFree_Error:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     : figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

[CHIRP,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select chirp file', root, 'MultiSelect','on');
CHIRP = cellstr(CHIRP)';

HeadFree = load(fullfile(root,CHIRP{1}),'GRAND','TRIAL','U','N');

filename = 'MakeFig_ChirpLog_HeadFree_Error';
xIdx = 1;

%% Frequency Domain
FIG = figure (1); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 3 2];
FIG.Name = filename;
movegui(FIG,'center')
hold on

errIdx = 4;
magIdx = 1;

amp = 3;
AMP  = HeadFree.U{1,3}{1}(amp);
FREQ = HeadFree.GRAND{amp,errIdx}.Mean{1}{7}(:,xIdx);
% VEL	= AMP*2*pi*FREQ;
ERROR = HeadFree.GRAND{amp,errIdx}.Mean{1}{8}(:,xIdx);
ESTD = HeadFree.GRAND{amp,errIdx}.STD{1}{8}(:,xIdx);

MAG = HeadFree.GRAND{amp,magIdx}.Mean{1}{8}(:,xIdx);
MSTD = 0*HeadFree.GRAND{amp,magIdx}.STD{1}{8}(:,xIdx);

ax1 = subplot(1,1,1);
    hold on
    ax1.Title.String = [num2str(AMP) , char(176)];
    ax1.Title.Color = 'w';
    ax1.Title.FontSize = 8;
    ax1.FontSize = 8;
    ax1.YLabel.String = ['Magnitude (' char(176) ')'];
    ax1.YLabel.FontSize = 8;
    ax1.YLim = [0 5.0];
    ax1.YTick = unique(sort([ax1.YTick ax1.YLim(2)]));
    ax1.XLabel.String = 'Frequency (Hz)';
    ax1.XLabel.FontSize = ax1.YLabel.FontSize;
    ax1.XLabel.Color = 'k';
    ax1.XLim = [0.1 12];
    ax1.XTick = [ax1.XLim(1),2:2:ax1.XLim(2)];
    
    h.patch = PlotPatch(MAG,   MSTD, FREQ, 0, HeadFree.N{1,1}, 'g', [0.4 0.4 0.6], 0.5, 2);

 	h.patch = PlotPatch(ERROR, ESTD, FREQ, 2, HeadFree.N{1,1}, 'y', [0.4 0.4 0.6], 0.5, 2);
 	
    ax1.XTick = unique(sort([min(ax1.XLim) ax1.XTick]));
    vel = round(AMP*2*pi*ax1.XTick);
    velLabel = cellfun(@(x) num2str(x), num2cell(vel), 'UniformOutput', false);

    ax3 = axes;
    ax3.Position = ax1.Position + [0 -0.00 0 0];
    ax3.FontSize = ax1.FontSize ;
    ax3.Color = 'none';
    ax3.YAxisLocation = 'right';
    ax3.YAxis.Color = 'none';
    ax3.XAxisLocation = 'top';
    ax3.XLim = ax1.XLim;
    ax3.XTick = ax1.XTick;
    ax3.XTickLabels = velLabel;
    ax3.XLabel.String = ['Peak Velocity (' char(176) '/s)'];
    ax3.XLabel.FontSize = ax1.YLabel.FontSize;

%% Time Domain
FIG = figure (2); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 8 2.5];
FIG.Name = filename;
movegui(FIG,'center')
hold on

errIdx = 4;
magIdx = 1;
amp = 3;

AMP  = HeadFree.U{1,3}{1}(amp);
TIME = HeadFree.GRAND{amp,errIdx}.Mean{1}{5}(:,xIdx);
ERROR = HeadFree.GRAND{amp,errIdx}.Mean{1}{6}(:,xIdx);
ESTD = HeadFree.GRAND{amp,errIdx}.STD{1}{6}(:,xIdx);
MAG = HeadFree.TRIAL{1,amp}{1,magIdx}.X(:,xIdx);
MSTD = 0*MAG;

ax1 = subplot(1,1,1);
    hold on
    ax1.Position(4) = ax1.Position(4)*0.8;
    ax1.Position = ax1.Position + [0 0.05 0 0];
    ax1.Title.String = [num2str(AMP) , char(176)];
    ax1.Title.Color = 'w';
    ax1.Title.FontSize = 8;
    ax1.FontSize = 8;
    ax1.YLabel.String = ['Magnitude (' char(176) ')'];
    ax1.YLabel.FontSize = 8;
    ax1.YLim = 20*[-1 1];
    ax1.YTick = unique(sort([ax1.YTick ax1.YLim(2)]));
    ax1.XLabel.String = 'Time (s)';
    ax1.XLabel.FontSize = ax1.YLabel.FontSize;
    ax1.XLabel.Color = 'k';
    ax1.XLim = [0 20];
    ax1.XTick = [ax1.XLim(1),2:2:ax1.XLim(2)];
    
    h.patch = PlotPatch(MAG,   MSTD, TIME, 0, HeadFree.N{1,1}, 'k', [0.4 0.4 0.6], 0.5, 2);

 	h.patch = PlotPatch(ERROR, ESTD, TIME, 2, HeadFree.N{1,1}, 'b', [0.4 0.4 0.6], 0.5, 2);
    
ax2 = axes;
    hold on
    ax2.FontSize = 8;
    ax2.Position = ax1.Position + [0 -0.0 0 0];
    ax2.Color = 'none';
    ax2.YColor = 'none';
    ax2.XLabel.String = 'Frequency (Hz)';
    ax2.XLabel.FontSize = 8;
    ax2.XAxisLocation = 'top';
    
    ax2.XLim = [0.1 12];
    ax2.XTick = [0.1 0.5 1 2 4 6 8 10 12];
    ax2.XScale = 'log';
    
%% Frequency Domain Log
FIG = figure (3); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 3 2];
FIG.Name = filename;
movegui(FIG,'center')
hold on

errIdx = 4;
magIdx = 1;

amp = 3;
AMP  = HeadFree.U{1,3}{1}(amp);
FREQ = HeadFree.GRAND{amp,errIdx}.Mean{1}{7}(:,xIdx);
ERROR = HeadFree.GRAND{amp,errIdx}.Mean{1}{8}(:,xIdx);
ESTD = HeadFree.GRAND{amp,errIdx}.STD{1}{8}(:,xIdx);

MAG = HeadFree.GRAND{amp,magIdx}.Mean{1}{8}(:,xIdx);
MSTD = 0*HeadFree.GRAND{amp,magIdx}.STD{1}{8}(:,xIdx);

ax1 = subplot(1,1,1);
    hold on
    ax1.Title.String = [num2str(AMP) , char(176)];
    ax1.Title.Color = 'w';
    ax1.Title.FontSize = 8;
    ax1.FontSize = 8;
    ax1.YLabel.String = ['Magnitude (' char(176) ')'];
    ax1.YLabel.FontSize = 8;
    ax1.YLim = [0 5.0];
    ax1.YTick = unique(sort([ax1.YTick ax1.YLim(2)]));
    ax1.XLabel.String = 'Frequency (Hz)';
    ax1.XLabel.FontSize = ax1.YLabel.FontSize;
    ax1.XLabel.Color = 'k';
    ax1.XLim = [0.1 12];
    ax1.XTick = [0.1 1 10];
  	ax1.XScale = 'log';
    ax1.XGrid = 'on';
    
    h.patch = PlotPatch(MAG,   MSTD, FREQ, 0, HeadFree.N{1,1}, 'g', [0.4 0.4 0.6], 0.5, 2);
    
 	h.patch = PlotPatch(ERROR(2:end), ESTD(2:end), FREQ(2:end), 2, HeadFree.N{1,1}, 'y', [0.4 0.4 0.6], 0.5, 2);
 	
    vel = round(AMP*2*pi*ax1.XTick);
    velLabel = cellfun(@(x) num2str(x), num2cell(vel), 'UniformOutput', false);

    ax3 = axes;
    ax3.Position = ax1.Position + [0 -0.00 0 0];
    ax3.FontSize = ax1.FontSize ;
    ax3.Color = 'none';
    ax3.YAxisLocation = 'right';
    ax3.YAxis.Color = 'none';
    ax3.XAxisLocation = 'top';
    ax3.XLim = ax1.XLim;
    ax3.XTick = ax1.XTick;
    ax3.XTickLabels = velLabel;
    ax3.XLabel.String = ['Peak Velocity (' char(176) '/s)'];
    ax3.XLabel.FontSize = ax1.YLabel.FontSize;
end