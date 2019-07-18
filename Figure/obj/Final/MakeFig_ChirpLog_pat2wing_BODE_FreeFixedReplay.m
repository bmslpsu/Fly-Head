function [FIG] = MakeFig_ChirpLog_pat2wing_BODE_FreeFixedReplay()
%% MakeFig_ChirpLog_pat2wing_BODE_FreeFixedReplay:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     : figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

[Free,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head-free file', root, 'MultiSelect','on');
Free = cellstr(Free)';

[Fixed,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head-fixed file', root, 'MultiSelect','on');
Fixed = cellstr(Fixed)';

[Replay,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select replay file', root, 'MultiSelect','on');
Replay = cellstr(Replay)';

HeadFree    = load(fullfile(root,Free{1}),'GRAND','U','N');
HeadFixed   = load(fullfile(root,Fixed{1}),'GRAND','U','N');
HeadReplay  = load(fullfile(root,Replay{1}),'GRAND','U','N');

%%
clearvars -except HeadFree HeadFixed HeadReplay

filename = 'ChirpLog_HeadFree_pat2head_BODE_ALL_new';
stim2wing = 8; 
err2wing = 6;
fixedhead = 3;
xIdx = 1;
CC = 'k';

% One Amplitude
FIG = figure (1); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 4 4];
FIG.Name = filename;
movegui(FIG,'center')
hold on

amp = 3;

HeadFree.AMP    = HeadFree.U{1,3}{1}(amp);
HeadFree.FREQ   = HeadFree.GRAND{amp,stim2wing}.Mean{2}{1}(:,xIdx);
HeadFree.VEL	= HeadFree.AMP*2*pi*HeadFree.FREQ;

HeadFree.GAIN  = HeadFree.GRAND{amp,stim2wing}.Mean{2}{2}(:,xIdx);
HeadFree.PHASE = rad2deg(HeadFree.GRAND{amp,stim2wing}.CircMean{7}{3}(:,xIdx));
HeadFree.GSTD  = HeadFree.GRAND{amp,stim2wing}.STD{2}{2}(:,xIdx);
HeadFree.PSTD  = rad2deg(HeadFree.GRAND{amp,stim2wing}.CircSTD{7}{3}(:,xIdx));

HeadFixed.GAIN  = HeadFixed.GRAND{1,fixedhead}.Mean{2}{2}(:,xIdx);
HeadFixed.PHASE = rad2deg(HeadFixed.GRAND{1,fixedhead}.CircMean{7}{3}(:,xIdx));
HeadFixed.GSTD  = HeadFixed.GRAND{1,fixedhead}.STD{2}{2}(:,xIdx);
HeadFixed.PSTD  = rad2deg(HeadFixed.GRAND{1,fixedhead}.CircSTD{7}{3}(:,xIdx));

HeadFree.OL_GAIN  = HeadFree.GRAND{amp,err2wing}.Mean{2}{2}(:,xIdx);
HeadFree.OL_PHASE = rad2deg(HeadFree.GRAND{amp,err2wing}.CircMean{7}{3}(:,xIdx));
HeadFree.OL_GSTD  = HeadFree.GRAND{amp,err2wing}.STD{2}{2}(:,xIdx);
HeadFree.OL_PSTD  = rad2deg(HeadFree.GRAND{amp,err2wing}.CircSTD{7}{3}(:,xIdx));

HeadReplay.GAIN  = HeadReplay.GRAND{1,fixedhead}.Mean{2}{2}(:,xIdx);
HeadReplay.PHASE = rad2deg(HeadReplay.GRAND{1,fixedhead}.CircMean{7}{3}(:,xIdx));
HeadReplay.GSTD  = HeadReplay.GRAND{1,fixedhead}.STD{2}{2}(:,xIdx);
HeadReplay.PSTD  = rad2deg(HeadReplay.GRAND{1,fixedhead}.CircSTD{7}{3}(:,xIdx));

[b,a] = butter(2,0.3,'low');
[bb,aa] = butter(2,0.3,'low');
mff = 1;

HeadFree.PHASE= filtfilt(bb,aa,medfilt1(filtfilt(b,a,HeadFree.PHASE),mff));
HeadFree.GAIN = filtfilt(bb,aa,medfilt1(filtfilt(b,a,HeadFree.GAIN),mff));
HeadFree.GSTD = filtfilt(bb,aa,medfilt1(filtfilt(b,a,HeadFree.GSTD),mff));
HeadFree.PSTD = filtfilt(bb,aa,medfilt1(filtfilt(b,a,HeadFree.PSTD),mff));

HeadFixed.PHASE= filtfilt(bb,aa,medfilt1(filtfilt(b,a,HeadFixed.PHASE),mff));
HeadFixed.GAIN = filtfilt(bb,aa,medfilt1(filtfilt(b,a,HeadFixed.GAIN),mff));
HeadFixed.GSTD = filtfilt(bb,aa,medfilt1(filtfilt(b,a,HeadFixed.GSTD),mff));
HeadFixed.PSTD = filtfilt(bb,aa,medfilt1(filtfilt(b,a,HeadFixed.PSTD),mff));

HeadFree.OL_PHASE= filtfilt(bb,aa,medfilt1(filtfilt(b,a,HeadFree.OL_PHASE),mff));
HeadFree.OL_GAIN = filtfilt(bb,aa,medfilt1(filtfilt(b,a,HeadFree.OL_GAIN),mff));
HeadFree.OL_GSTD = filtfilt(bb,aa,medfilt1(filtfilt(b,a,HeadFree.OL_GSTD),mff));
HeadFree.OL_PSTD = filtfilt(bb,aa,medfilt1(filtfilt(b,a,HeadFree.OL_PSTD),mff));

HeadReplay.PHASE= filtfilt(bb,aa,medfilt1(filtfilt(b,a,HeadReplay.PHASE),mff));
HeadReplay.GAIN = filtfilt(bb,aa,medfilt1(filtfilt(b,a,HeadReplay.GAIN),mff));
HeadReplay.GSTD = filtfilt(bb,aa,medfilt1(filtfilt(b,a,HeadReplay.GSTD),mff));
HeadReplay.PSTD = filtfilt(bb,aa,medfilt1(filtfilt(b,a,HeadReplay.PSTD),mff));

ax1 = subplot(2,1,1);
    hold on
    ax1.Title.String = [num2str(HeadFree.AMP) , char(176)];
    ax1.Title.Color = 'w';
    ax1.Title.FontSize = 8;
    ax1.FontSize = 8;
    ax1.YLabel.String = ['Gain (' char(176) '/' char(176) ')'];
    ax1.YLabel.FontSize = 8;
    ax1.YLim = [0 0.5];
    ax1.YTick = unique(sort([ax1.YTick ax1.YLim(2)]));
    ax1.XLabel.String = 'Frequency (Hz)';
    ax1.XLabel.FontSize = ax1.YLabel.FontSize;
    ax1.XLabel.Color = 'w';
    ax1.XLim = [0.3 10];
    ax1.XTick = [ax1.XLim(1),2:2:ax1.XLim(2)];
    ax1.XTickLabels = '';
    
    [~,h.stim2wing] = PlotPatch(HeadFree.GAIN, HeadFree.GSTD, HeadFree.FREQ,...
                        2,HeadFree.N{1,1},'r',[0.4 0.4 0.6],0.5,2);
                    
    [~,h.err2wing] = PlotPatch(HeadFree.OL_GAIN, HeadFree.OL_GSTD, HeadFree.FREQ,...
                        2,HeadFree.N{1,1},'c',[0.4 0.4 0.6],0.5,2);
                    
    [~,h.fixed] = PlotPatch(HeadFixed.GAIN, HeadFixed.GSTD, HeadFree.FREQ,...
                        2,HeadFixed.N{1,1},[0.4 0 0.8],[0.4 0.4 0.6],0.5,2);
	
    [~,h.replay] = PlotPatch(HeadReplay.GAIN, HeadReplay.GSTD, HeadFree.FREQ,...
                        2,HeadReplay.N{1,1},'y',[0.4 0.4 0.6],0.5,2);
	
  	uistack([h.stim2wing h.err2wing],'top')

    ax1.XTick = unique(sort([min(ax1.XLim) ax1.XTick]));
    vel = round(HeadFree.AMP*2*pi*ax1.XTick);
    velLabel = cellfun(@(x) num2str(x), num2cell(vel), 'UniformOutput', false);

ax2 = subplot(2,1,2);
    hold on
%     ax2.Title.String = [num2str(HeadFree.AMP) , char(176)];
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
    
    [~,h.stim2wing] = PlotPatch(HeadFree.PHASE, HeadFree.PSTD,HeadFree.FREQ,...
                        2, HeadFree.N{1,1}, 'r', [0.4 0.4 0.6], 0.5, 2);
  	
    [~,h.err2wing] = PlotPatch(HeadFree.OL_PHASE, HeadFree.OL_PSTD,HeadFree.FREQ,...
                        2, HeadFree.N{1,1}, 'c', [0.4 0.4 0.6], 0.5, 2);
 	
    [~,h.fixed] = PlotPatch(HeadFixed.PHASE, HeadFixed.PSTD,HeadFree.FREQ,...
                        2, HeadFixed.N{1,1}, [0.4 0 0.8], [0.4 0.4 0.6], 0.5, 2);
    
    [~,h.replay] = PlotPatch(HeadReplay.PHASE, HeadReplay.PSTD,HeadFree.FREQ,...
                        2, HeadReplay.N{1,1}, 'y', [0.4 0.4 0.6], 0.5, 2);
                    
    uistack([h.stim2wing h.err2wing],'top')
    
    plot([0 12],[0 0],'--k','LineWidth',1);

    ax2.XTick = ax1.XTick;

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