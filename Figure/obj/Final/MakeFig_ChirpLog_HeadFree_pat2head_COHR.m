function [FIG] = MakeFig_ChirpLog_HeadFree_pat2head_COHR()
%% MakeFig_ChirpLog_HeadFree_pat2head_COHR_all:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     : figure handle
%

root = 'H:\DATA\Rigid_Data\';

[Free,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head-free chirp file', root, 'MultiSelect','on');
Free = cellstr(Free)';

[Fixed,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head-fixed chirp file', root, 'MultiSelect','on');
Fixed = cellstr(Fixed)';

HeadFree = load(fullfile(root,Free{1}),'GRAND','U','N');
HeadFixed = load(fullfile(root,Fixed{1}),'GRAND','U','N');

%%
figNum = 1;
filename = 'ChirpLog_HeadFree_pat2head_COHR_all';
xIdx = 1;
CC = [0 0 1];
offset = 0.05;

FIG = figure (figNum); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 4 3];
FIG.Name = filename;
movegui(FIG,'center')
hold on

amp = 1;
Amp = HeadFree.U{1,3}{1}(amp);
ax1 = gca;
    hold on
    ax1.Position = ax1.Position + [0 offset/2 0 -offset];
% 	ax1.Title.String = [num2str(Amp) , char(176)];
    ax1.Title.Color = 'k';
    ax1.Title.FontSize = 8;
    ax1.FontSize = 8;
    ax1.YLabel.String = 'Coherence';
    ax1.YLabel.FontSize = 8;
    ax1.YLim = [0 1.0];
    ax1.YTick = unique(sort([ax1.YTick ax1.YLim(2)]));
    ax1.XLabel.String = 'Frequency (Hz)';
    ax1.XLabel.Color = 'k';
    ax1.XLabel.FontSize = ax1.YLabel.FontSize;
    ax1.XLim = [0.1 12];

    FREQ = HeadFree.GRAND{amp,5}.Mean{1}{8};
    COHR = HeadFree.GRAND{amp,5}.Mean{1}{7}(:,xIdx);
    STD  = HeadFree.GRAND{amp,5}.STD{1}{7}(:,xIdx);

    [~,h.head] = PlotPatch(COHR, STD, FREQ , 1, HeadFree.N{1,1}, CC, [0.4 0.4 0.6], 0.5, 1);
    
	FREQ = HeadFree.GRAND{amp,8}.Mean{1}{8};
    COHR = HeadFree.GRAND{amp,8}.Mean{1}{7}(:,xIdx);
    STD  = HeadFree.GRAND{amp,8}.STD{1}{7}(:,xIdx);

    [~,h.wing] = PlotPatch(COHR, STD, FREQ , 1, HeadFree.N{1,1}, 'r', [0.4 0.4 0.6], 0.5, 1);
    
  	FREQ = HeadFixed.GRAND{amp,3}.Mean{1}{8};
    COHR = HeadFixed.GRAND{amp,3}.Mean{1}{7}(:,xIdx);
    STD  = HeadFixed.GRAND{amp,3}.STD{1}{7}(:,xIdx);

    [~,h.wing_fixed] = PlotPatch(COHR, STD, FREQ , 1, HeadFree.N{1,1}, [0.3 0.1 0.7], [0.4 0.4 0.6], 0.5, 1);
    
    uistack(h.head,'top')
    uistack(h.wing,'top')
    
    ax1.XTick = unique(sort([min(ax1.XLim) ax1.XTick]));
    vel = round(Amp*2*pi*ax1.XTick);
    velLabel = cellfun(@(x) num2str(x), num2cell(vel), 'UniformOutput', false);               
    
    ax3 = axes;
    ax3.Position = ax1.Position;
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