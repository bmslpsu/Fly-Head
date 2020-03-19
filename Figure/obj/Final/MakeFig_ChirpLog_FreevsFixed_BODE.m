function [FIG] = MakeFig_ChirpLog_FreevsFixed_BODE()
%% MakeFig_ChirpLog_FreevsFixed_BODE:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     : figure handle
%

root = 'H:\DATA\Rigid_Data\';

[CHIRP,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select chirp file', root, 'MultiSelect','on');
CHIRP = cellstr(CHIRP)';

HeadFree = load(fullfile(root,CHIRP{1}),'GRAND','U','N');
clearvars -except HeadFree

%%
filename = 'ChirpLog_FreevsFixed_BODE';
catIdx = 6;
xIdx = 1;
y_lim = 0.3;
n_amp = HeadFree.N.Amp;

AMP     = [];
FREQ    = [];
VEL     = [];
GAIN    = [];
PHASE   = [];
for jj = 1:n_amp % amplitudes
	AMP(:,jj)      = HeadFree.U{1,3}{1}(jj);
    FREQ(:,jj)     = HeadFree.GRAND{jj,catIdx}.Mean{2}{1}(:,xIdx);
    VEL(:,jj)      = AMP(:,jj)*2*pi*FREQ(:,jj);
    GAIN(:,jj)     = abs(HeadFree.GRAND{jj,catIdx}.Mean{1}{17}(:,xIdx));
    PHASE(:,jj)    = rad2deg(angle(HeadFree.GRAND{jj,catIdx}.Mean{1}{17}(:,xIdx)));
    GSTD(:,jj)     = abs(HeadFree.GRAND{jj,catIdx}.STD{1}{17}(:,xIdx));
    PSTD(:,jj)     = rad2deg((HeadFree.GRAND{jj,catIdx}.STD{1}{17}(:,xIdx)));
    
 	[b,a] = butter(2,0.4,'low');
    [bb,aa] = butter(2,0.4,'low');
    mff = 2;
    PHASE(:,jj) = filtfilt(bb,aa,medfilt1(filtfilt(b,a,PHASE(:,jj)),mff));
    GAIN(:,jj)  = filtfilt(bb,aa,medfilt1(filtfilt(b,a,GAIN(:,jj)),mff));
    GSTD(:,jj)  = filtfilt(bb,aa,medfilt1(filtfilt(b,a,GSTD(:,jj)),mff));
    PSTD(:,jj)  = filtfilt(bb,aa,medfilt1(filtfilt(b,a,PSTD(:,jj)),mff));
end
flip = (PHASE<-60) & (FREQ<0.4);
PHASE(flip) = PHASE(flip) + 360;

%% 15 deg Amplitude
FIG = figure (1); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 1.5*2.5 5];
FIG.Name = filename;
movegui(FIG,'center')
CC = hsv(n_amp+1);
pp = 1;
ax = gobjects(3,1);
amp = 3;
for jj = amp % amplitudes    
    ax(1,pp) = subplot(2,1,pp);
        hold on
        ax(1,pp).YLabel.String = 'Gain (°/°)';
        ax(1,pp).YLabel.FontSize = 10;
        ax(1,pp).YLim = [0 y_lim];
        ax(1,pp).YTick = unique(sort([ax(1,pp).YTick ax(1,pp).YLim(2)]));

        h.patch = PlotPatch(GAIN(:,jj), GSTD(:,jj), FREQ(:,jj) ,...
                            1, HeadFree.N{1,1}, CC(jj+1,:), [0.4 0.4 0.6], 0.5, 1);
                               
        ax(1,pp).XTick = unique(sort([min(ax(1,pp).XLim) 2:2:12]));
        vel = round(AMP(:,jj)*2*pi*ax(1,pp).XTick);
        velLabel = cellfun(@(x) num2str(x), num2cell(vel), 'UniformOutput', false);
               
    ax(2,pp) = subplot(2,1,pp + 1);
        hold on
        ax(2,pp).Title.String = [num2str(AMP(:,jj)) , char(176)];
        ax(2,pp).Title.Color = 'k';
        ax(2,pp).Title.FontSize = ax(1,pp).Title.FontSize;
        ax(2,pp).YLabel.String = 'Phase Difference (°)';
        ax(2,pp).YLabel.FontSize = ax(1,pp).YLabel.FontSize;
        ax(2,pp).XLabel.String = 'Frequency (Hz)';
        ax(2,pp).XLabel.FontSize = ax(1,pp).YLabel.FontSize;
      	ax(2,pp).XLabel.Color = 'k';
        ax(2,pp).YLim = rad2deg(pi*[-1 1]);
        ax(2,pp).YTick = -180:60:180;
        
        h.patch = PlotPatch(PHASE(:,jj), PSTD(:,jj),FREQ(:,jj),...
            1, HeadFree.N{1,1} ,CC(jj+1,:), [0.4 0.4 0.6],0.5,1);
                
        plot([0 12],[0 0],'--k','LineWidth',1)
                
	ax(3,pp) = axes;
        ax(3,pp).Position = ax(1,pp).Position + [0 -0.00 0 0];
        ax(3,pp).FontSize = ax(1,pp).FontSize ;
        ax(3,pp).Color = 'none';
        ax(3,pp).YAxisLocation = 'right';
        ax(3,pp).YAxis.Color = 'none';
        ax(3,pp).XAxisLocation = 'top';
        ax(3,pp).XLim = ax(1,pp).XLim;
        ax(3,pp).XTick = ax(1,pp).XTick;
        ax(3,pp).XTickLabels = velLabel;
      	ax(3,pp).XLabel.String = 'Peak Velocity (°/s)';
        ax(3,pp).XLabel.FontSize = ax(1,pp).YLabel.FontSize;
        
  	pp = pp + 1;
end
set(ax,'FontSize', 8, 'LineWidth', 1.5, 'XLim',[0.1 12])
set(ax, 'XTick', sort([ax(1,1).XLim(1),2:2:12]))
linkaxes(ax,'x')

end