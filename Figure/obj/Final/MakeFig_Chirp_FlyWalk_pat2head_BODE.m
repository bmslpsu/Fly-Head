function [FIG] = MakeFig_Chirp_FlyWalk_pat2head_BODE()
%% MakeFig_Chirp_FlyWalk_pat2head_BODE:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     : figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

[flying,path1] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select chirp file', root, 'MultiSelect','on');
flying = cellstr(flying)';

[walking,path2] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select chirp file', root, 'MultiSelect','on');
walking = cellstr(walking)';

FLY = load(fullfile(path1,flying{1}),'GRAND','U','N');
WALK = load(fullfile(path2,walking{1}),'GRAND','U','N');

%%
figNum = 1;
filename = 'Chirp_FlyWalk_pat2head_BODE';
catIdx = 5;
xIdx = 1;
CC = [0 0 0.7];

FIG = figure (figNum); clf
FIG.Color = 'w';
FIG.Position = [100 100 650 750];
FIG.Name = filename;
movegui(FIG,'center')
hold on

[b,a] = butter(2,0.3,'low');
mI = 3;
for jj = 1:FLY.N{1,3} % amplitudes
	AMP(:,jj)   = FLY.U{1,3}{1}(jj);
    FREQ(:,jj)  = FLY.GRAND{jj,catIdx}.Mean{2}{1}(:,xIdx);
    VEL(:,jj) 	= AMP(:,jj)*2*pi*FREQ(:,jj);
    
    GAIN(:,jj)     = filtfilt(b,a,medfilt1(FLY.GRAND{jj,catIdx}.Mean{2}{2}(:,xIdx),mI));
    PHASE(:,jj)    = filtfilt(b,a,medfilt1(rad2deg(FLY.GRAND{jj,catIdx}.CircMean{7}{3}(:,xIdx)),mI));
 	STD(:,jj)      = filtfilt(b,a,medfilt1(FLY.GRAND{jj,catIdx}.STD{2}{2}(:,xIdx),mI));
    PSTD(:,jj)     = filtfilt(b,a,medfilt1(rad2deg(FLY.GRAND{jj,catIdx}.CircSTD{7}{3}(:,xIdx)),mI));   
end

fly.FREQ    = mean(FREQ,2);
fly.GAIN    = mean(GAIN,2);
fly.STD     = std(GAIN,[],2);
fly.PHASE  	= mean(PHASE,2);
fly.PSTD   	= std(PHASE,[],2);

walk.FREQ 	= WALK.GRAND{1,4}.Mean{2}{1}(:,xIdx);
walk.GAIN   = filtfilt(b,a,medfilt1(WALK.GRAND{1,4}.Mean{2}{2}(:,xIdx),mI));
walk.STD    = filtfilt(b,a,medfilt1(WALK.GRAND{1,4}.STD{2}{2}(:,xIdx),mI));
walk.PHASE	= rad2deg(filtfilt(b,a,medfilt1(WALK.GRAND{1,4}.CircMean{7}{3}(:,xIdx),mI)));
walk.PSTD  	= rad2deg(filtfilt(b,a,medfilt1(WALK.GRAND{1,4}.CircSTD{7}{3}(:,xIdx),mI)));

ax1 = subplot(2,1,1);
    hold on
    ax1.Title.String = [num2str(AMP(:,jj)) , char(176)];
    ax1.Title.Color = 'w';
    ax1.Title.FontSize = 16;
    ax1.FontSize = 12;
    ax1.YLabel.String = ['Gain (' char(176) '/' char(176) ')'];
    ax1.YLabel.FontSize = 14;
    ax1.YLim = [0 1.2];
    ax1.YTick = unique(sort([ax1.YTick ax1.YLim(2)]));
    ax1.XLim = [0.3 10];
 	ax1.XTickLabels = '';

    [~,h.fly] = PlotPatch(fly.GAIN, fly.STD, fly.FREQ  , 2, FLY.N{1,1}, CC, [0.4 0.4 0.6], 0.5, 3);

    [~,h.walk] = PlotPatch(walk.GAIN, walk.STD, walk.FREQ , 2, WALK.N{1,1}, [0.5 0 0], [0.4 0.4 0.6], 0.5, 3);

    Leg = legend([h.fly h.walk], 'Flying', 'Walking');
    Leg.FontSize = 14;

    ax1.XTick = unique(sort([min(ax1.XLim) ax1.XTick]));

ax2 = subplot(2,1,2);
    hold on
    ax2.Title.Color = 'k';
    ax2.Title.FontSize = ax1.Title.FontSize;
    ax2.FontSize = ax1.FontSize;
    ax2.YLabel.String = ['Phase Difference (' char(176) ')'];
    ax2.YLabel.FontSize = ax1.YLabel.FontSize;
    ax2.XLabel.String = 'Frequency (Hz)';
    ax2.XLabel.FontSize = ax1.YLabel.FontSize;
    ax2.XLabel.Color = 'k';
    ax2.YLim = rad2deg(pi*[-1 1]);
    ax2.XLim = ax1.XLim;
    ax2.YTick = -180:60:180;

    [~,h.fly] = PlotPatch(fly.PHASE, PSTD(:,2), fly.FREQ  , 2, FLY.N{1,1}, CC, [0.4 0.4 0.6], 0.5, 3);

    [~,h.walk] = PlotPatch(walk.PHASE, walk.PSTD, walk.FREQ , 2, WALK.N{1,1}, [0.5 0 0], [0.4 0.4 0.6], 0.5, 3);

    plot([0 12],[0 0],'--g','LineWidth',1);

    ax2.XTick = ax1.XTick;

end