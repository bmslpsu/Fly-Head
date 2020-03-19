function [FIG] = MakeFig_ChirpLog_HeadFreeFixed_wing_BODE()
%% MakeFig_ChirpLog_HeadFreeFixed_wing_BODE:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     : figure handle
%

root = 'H:\DATA\Rigid_Data\';

[Free,~] = uigetfile({'*.mat'}, 'Select free file', root, 'MultiSelect','off');
[Fixed,~] = uigetfile({'*.mat'}, 'Select fixed file', root, 'MultiSelect','off');

HeadFree = load(fullfile(root,Free),'GRAND','U','N');
HeadFixed = load(fullfile(root,Fixed),'GRAND','U','N');

clearvars -except HeadFree

%%
filename = 'ChirpLog_HeadFree_pat2wing_BODE_ALL_new';
pat2wing_free = 8;
err2wing_free = 6;
pat2wing_fixed = 3;

xIdx = 1;
y_lim = 0.5;
n_amp = HeadFree.N.Amp;

%% All Amplitudes
FIG = figure (1); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 n_amp*2.5 5];
FIG.Name = filename;
movegui(FIG,'center')
CC = hsv(n_amp+1);

%% pat2wing , free
AMP     = [];
FREQ    = [];
VEL     = [];
GAIN    = [];
PHASE   = [];
for jj = 1:n_amp % amplitudes
	AMP(:,jj)      = HeadFree.U{1,3}{1}(jj);
    FREQ(:,jj)     = HeadFree.GRAND{jj,pat2wing_free}.Mean{2}{1}(:,xIdx);
    VEL(:,jj)      = AMP(:,jj)*2*pi*FREQ(:,jj);

	% GAIN(:,jj)     = abs(HeadFree.GRAND{jj,pat2wing_free}.Mean{1}{17}(:,xIdx));
    % PHASE(:,jj)    = rad2deg(angle(HeadFree.GRAND{jj,pat2wing_free}.Mean{1}{17}(:,xIdx)));
 	GAIN(:,jj)     = HeadFree.GRAND{jj,pat2wing_free}.Mean{1}{2}(:,xIdx);
    PHASE(:,jj)    = rad2deg(HeadFree.GRAND{jj,pat2wing_free}.CircMean{7}{3}(:,xIdx));
    
    % GSTD(:,jj)     = abs(HeadFree.GRAND{jj,pat2wing_free}.STD{1}{17}(:,xIdx));
    % PSTD(:,jj)     = rad2deg((HeadFree.GRAND{jj,pat2wing_free}.STD{1}{17}(:,xIdx)));
    
  	GSTD(:,jj)     = HeadFree.GRAND{jj,pat2wing_free}.STD{1}{2}(:,xIdx);
    PSTD(:,jj)     = rad2deg(HeadFree.GRAND{jj,pat2wing_free}.CircSTD{7}{3}(:,xIdx));
    
 	[b,a] = butter(2,0.5,'low');
    [bb,aa] = butter(2,0.5,'low');
    mff = 2;
    PHASE(:,jj) = filtfilt(bb,aa,medfilt1(filtfilt(b,a,PHASE(:,jj)),mff));
    GAIN(:,jj)  = filtfilt(bb,aa,medfilt1(filtfilt(b,a,GAIN(:,jj)),mff));
    GSTD(:,jj)  = filtfilt(bb,aa,medfilt1(filtfilt(b,a,GSTD(:,jj)),mff));
    PSTD(:,jj)  = filtfilt(bb,aa,medfilt1(filtfilt(b,a,PSTD(:,jj)),mff));
end
% flip = (PHASE<-60) & (FREQ<0.4);
% PHASE(flip) = PHASE(flip) + 360;

pp = 1;
ax = gobjects(3,n_amp);
for jj = 1:n_amp % amplitudes    
    ax(1,jj) = subplot(2,n_amp,pp);
        hold on
        ax(1,jj).YLabel.String = 'Gain (°/°)';
        ax(1,jj).YLabel.FontSize = 10;
        ax(1,jj).YLim = [0 y_lim];
        ax(1,jj).YTick = unique(sort([ax(1,jj).YTick ax(1,jj).YLim(2)]));

        h.patch = PlotPatch(GAIN(:,jj), GSTD(:,jj), FREQ(:,jj) ,...
                            1, HeadFree.N{1,1}, 'g', [0.4 0.4 0.6], 0.5, 1);
                               
        ax(1,jj).XTick = unique(sort([min(ax(1,jj).XLim) 2:2:12]));
        vel = round(AMP(:,jj)*2*pi*ax(1,jj).XTick);
        velLabel = cellfun(@(x) num2str(x), num2cell(vel), 'UniformOutput', false);
               
    ax(2,jj) = subplot(2,HeadFree.N{1,3},pp + HeadFree.N{1,3});
        hold on
        ax(2,jj).Title.String = [num2str(AMP(:,jj)) , char(176)];
        ax(2,jj).Title.Color = 'k';
        ax(2,jj).Title.FontSize = ax(1,jj).Title.FontSize;
        ax(2,jj).YLabel.String = 'Phase Difference (°)';
        ax(2,jj).YLabel.FontSize = ax(1,jj).YLabel.FontSize;
        ax(2,jj).XLabel.String = 'Frequency (Hz)';
        ax(2,jj).XLabel.FontSize = ax(1,jj).YLabel.FontSize;
      	ax(2,jj).XLabel.Color = 'k';
        ax(2,jj).YLim = rad2deg(pi*[-1 1]);
        ax(2,jj).YTick = -180:60:180;
        
        h.patch = PlotPatch(PHASE(:,jj), PSTD(:,jj),FREQ(:,jj),...
            1, HeadFree.N{1,1} , 'g', [0.4 0.4 0.6],0.5,1);
                
        plot([0 12],[0 0],'--k','LineWidth',1)
                
	ax(3,jj) = axes;
        ax(3,jj).Position = ax(1,jj).Position + [0 -0.00 0 0];
        ax(3,jj).FontSize = ax(1,jj).FontSize ;
        ax(3,jj).Color = 'none';
        ax(3,jj).YAxisLocation = 'right';
        ax(3,jj).YAxis.Color = 'none';
        ax(3,jj).XAxisLocation = 'top';
        ax(3,jj).XLim = ax(1,jj).XLim;
        ax(3,jj).XTick = ax(1,jj).XTick;
        ax(3,jj).XTickLabels = velLabel;
      	ax(3,jj).XLabel.String = 'Peak Velocity (°/s)';
        ax(3,jj).XLabel.FontSize = ax(1,jj).YLabel.FontSize;
        
  	pp = pp + 1;
end
set(ax,'FontSize', 8, 'LineWidth', 1.5, 'XLim',[0.1 12])
set(ax, 'XTick', sort([ax(1,jj).XLim(1),2:2:12]))
linkaxes(ax,'x')


%% err2wing , free
AMP     = [];
FREQ    = [];
VEL     = [];
GAIN    = [];
PHASE   = [];
for jj = 1:n_amp % amplitudes
	AMP(:,jj)      = HeadFree.U{1,3}{1}(jj);
    FREQ(:,jj)     = HeadFree.GRAND{jj,err2wing_free}.Mean{2}{1}(:,xIdx);
    VEL(:,jj)      = AMP(:,jj)*2*pi*FREQ(:,jj);

	% GAIN(:,jj)     = abs(HeadFree.GRAND{jj,err2wing_free}.Mean{1}{17}(:,xIdx));
    % PHASE(:,jj)    = rad2deg(angle(HeadFree.GRAND{jj,err2wing_free}.Mean{1}{17}(:,xIdx)));
 	GAIN(:,jj)     = HeadFree.GRAND{jj,err2wing_free}.Mean{1}{2}(:,xIdx);
    PHASE(:,jj)    = rad2deg(HeadFree.GRAND{jj,err2wing_free}.CircMean{7}{3}(:,xIdx));
    
    % GSTD(:,jj)     = abs(HeadFree.GRAND{jj,err2wing_free}.STD{1}{17}(:,xIdx));
    % PSTD(:,jj)     = rad2deg((HeadFree.GRAND{jj,err2wing_free}.STD{1}{17}(:,xIdx)));
    
  	GSTD(:,jj)     = HeadFree.GRAND{jj,err2wing_free}.STD{1}{2}(:,xIdx);
    PSTD(:,jj)     = rad2deg(HeadFree.GRAND{jj,err2wing_free}.CircSTD{7}{3}(:,xIdx));
    
 	[b,a] = butter(2,0.5,'low');
    [bb,aa] = butter(2,0.5,'low');
    mff = 2;
    PHASE(:,jj) = filtfilt(bb,aa,medfilt1(filtfilt(b,a,PHASE(:,jj)),mff));
    GAIN(:,jj)  = filtfilt(bb,aa,medfilt1(filtfilt(b,a,GAIN(:,jj)),mff));
    GSTD(:,jj)  = filtfilt(bb,aa,medfilt1(filtfilt(b,a,GSTD(:,jj)),mff));
    PSTD(:,jj)  = filtfilt(bb,aa,medfilt1(filtfilt(b,a,PSTD(:,jj)),mff));
end
% flip = (PHASE<-60) & (FREQ<0.4);
% PHASE(flip) = PHASE(flip) + 360;

pp = 1;
ax = gobjects(3,n_amp);
for jj = 1:n_amp % amplitudes    
    ax(1,jj) = subplot(2,n_amp,pp);
        hold on
        ax(1,jj).YLabel.String = 'Gain (°/°)';
        ax(1,jj).YLabel.FontSize = 10;
        ax(1,jj).YLim = [0 y_lim];
        ax(1,jj).YTick = unique(sort([ax(1,jj).YTick ax(1,jj).YLim(2)]));

        h.patch = PlotPatch(GAIN(:,jj), GSTD(:,jj), FREQ(:,jj) ,...
                            1, HeadFree.N{1,1}, 'c', [0.4 0.4 0.6], 0.5, 1);
                               
        ax(1,jj).XTick = unique(sort([min(ax(1,jj).XLim) 2:2:12]));
        vel = round(AMP(:,jj)*2*pi*ax(1,jj).XTick);
        velLabel = cellfun(@(x) num2str(x), num2cell(vel), 'UniformOutput', false);
               
    ax(2,jj) = subplot(2,HeadFree.N{1,3},pp + HeadFree.N{1,3});
        hold on
        ax(2,jj).Title.String = [num2str(AMP(:,jj)) , char(176)];
        ax(2,jj).Title.Color = 'k';
        ax(2,jj).Title.FontSize = ax(1,jj).Title.FontSize;
        ax(2,jj).YLabel.String = 'Phase Difference (°)';
        ax(2,jj).YLabel.FontSize = ax(1,jj).YLabel.FontSize;
        ax(2,jj).XLabel.String = 'Frequency (Hz)';
        ax(2,jj).XLabel.FontSize = ax(1,jj).YLabel.FontSize;
      	ax(2,jj).XLabel.Color = 'k';
        ax(2,jj).YLim = rad2deg(pi*[-1 1]);
        ax(2,jj).YTick = -180:60:180;
        
        h.patch = PlotPatch(PHASE(:,jj), PSTD(:,jj),FREQ(:,jj),...
            1, HeadFree.N{1,1} ,'c', [0.4 0.4 0.6],0.5,1);
                
        plot([0 12],[0 0],'--k','LineWidth',1)
                
	ax(3,jj) = axes;
        ax(3,jj).Position = ax(1,jj).Position + [0 -0.00 0 0];
        ax(3,jj).FontSize = ax(1,jj).FontSize ;
        ax(3,jj).Color = 'none';
        ax(3,jj).YAxisLocation = 'right';
        ax(3,jj).YAxis.Color = 'none';
        ax(3,jj).XAxisLocation = 'top';
        ax(3,jj).XLim = ax(1,jj).XLim;
        ax(3,jj).XTick = ax(1,jj).XTick;
        ax(3,jj).XTickLabels = velLabel;
      	ax(3,jj).XLabel.String = 'Peak Velocity (°/s)';
        ax(3,jj).XLabel.FontSize = ax(1,jj).YLabel.FontSize;
        
  	pp = pp + 1;
end
set(ax,'FontSize', 8, 'LineWidth', 1.5, 'XLim',[0.1 12])
set(ax, 'XTick', sort([ax(1,jj).XLim(1),2:2:12]))
linkaxes(ax,'x')


%% err2wing , fixed
AMP     = [];
FREQ    = [];
VEL     = [];
GAIN    = [];
PHASE   = [];
for jj = 1:n_amp % amplitudes
	AMP(:,jj)      = HeadFixed.U{1,3}{1}(jj);
    FREQ(:,jj)     = HeadFixed.GRAND{jj,pat2wing_fixed}.Mean{2}{1}(:,xIdx);
    VEL(:,jj)      = AMP(:,jj)*2*pi*FREQ(:,jj);

	% GAIN(:,jj)     = abs(HeadFixed.GRAND{jj,pat2wing_fixed}.Mean{1}{17}(:,xIdx));
    % PHASE(:,jj)    = rad2deg(angle(HeadFixed.GRAND{jj,pat2wing_fixed}.Mean{1}{17}(:,xIdx)));
 	GAIN(:,jj)     = HeadFixed.GRAND{jj,pat2wing_fixed}.Mean{1}{2}(:,xIdx);
    PHASE(:,jj)    = rad2deg(HeadFixed.GRAND{jj,pat2wing_fixed}.CircMean{7}{3}(:,xIdx));
    
    % GSTD(:,jj)     = abs(HeadFixed.GRAND{jj,pat2wing_fixed}.STD{1}{17}(:,xIdx));
    % PSTD(:,jj)     = rad2deg((HeadFixed.GRAND{jj,pat2wing_fixed}.STD{1}{17}(:,xIdx)));
    
  	GSTD(:,jj)     = HeadFixed.GRAND{jj,pat2wing_fixed}.STD{1}{2}(:,xIdx);
    PSTD(:,jj)     = rad2deg(HeadFixed.GRAND{jj,pat2wing_fixed}.CircSTD{7}{3}(:,xIdx));
    
 	[b,a] = butter(2,0.5,'low');
    [bb,aa] = butter(2,0.5,'low');
    mff = 2;
    PHASE(:,jj) = filtfilt(bb,aa,medfilt1(filtfilt(b,a,PHASE(:,jj)),mff));
    GAIN(:,jj)  = filtfilt(bb,aa,medfilt1(filtfilt(b,a,GAIN(:,jj)),mff));
    GSTD(:,jj)  = filtfilt(bb,aa,medfilt1(filtfilt(b,a,GSTD(:,jj)),mff));
    PSTD(:,jj)  = filtfilt(bb,aa,medfilt1(filtfilt(b,a,PSTD(:,jj)),mff));
    endTD(:,jj)  = filtfilt(bb,aa,medfilt1(filtfilt(b,a,PSTD(:,jj)),mff));
end
% flip = (PHASE<-60) & (FREQ<0.4);
% PHASE(flip) = PHASE(flip) + 360;

pp = 1;
ax = gobjects(3,n_amp);
for jj = 1:n_amp % amplitudes    
    ax(1,jj) = subplot(2,n_amp,pp);
        hold on
        ax(1,jj).YLabel.String = 'Gain (°/°)';
        ax(1,jj).YLabel.FontSize = 10;
        ax(1,jj).YLim = [0 y_lim];
        ax(1,jj).YTick = unique(sort([ax(1,jj).YTick ax(1,jj).YLim(2)]));

        h.patch = PlotPatch(GAIN(:,jj), GSTD(:,jj), FREQ(:,jj) ,...
                            1, HeadFixed.N{1,1}, 'r', [0.4 0.4 0.6], 0.5, 1);
                               
        ax(1,jj).XTick = unique(sort([min(ax(1,jj).XLim) 2:2:12]));
        vel = round(AMP(:,jj)*2*pi*ax(1,jj).XTick);
        velLabel = cellfun(@(x) num2str(x), num2cell(vel), 'UniformOutput', false);
               
    ax(2,jj) = subplot(2,HeadFixed.N{1,3},pp + HeadFixed.N{1,3});
        hold on
        ax(2,jj).Title.String = [num2str(AMP(:,jj)) , char(176)];
        ax(2,jj).Title.Color = 'k';
        ax(2,jj).Title.FontSize = ax(1,jj).Title.FontSize;
        ax(2,jj).YLabel.String = 'Phase Difference (°)';
        ax(2,jj).YLabel.FontSize = ax(1,jj).YLabel.FontSize;
        ax(2,jj).XLabel.String = 'Frequency (Hz)';
        ax(2,jj).XLabel.FontSize = ax(1,jj).YLabel.FontSize;
      	ax(2,jj).XLabel.Color = 'k';
        ax(2,jj).YLim = rad2deg(pi*[-1 1]);
        ax(2,jj).YTick = -180:60:180;
        
        h.patch = PlotPatch(PHASE(:,jj), PSTD(:,jj),FREQ(:,jj),...
            1, HeadFixed.N{1,1} ,'r', [0.4 0.4 0.6],0.5,1);
                
        plot([0 12],[0 0],'--k','LineWidth',1)
                
	ax(3,jj) = axes;
        ax(3,jj).Position = ax(1,jj).Position + [0 -0.00 0 0];
        ax(3,jj).FontSize = ax(1,jj).FontSize ;
        ax(3,jj).Color = 'none';
        ax(3,jj).YAxisLocation = 'right';
        ax(3,jj).YAxis.Color = 'none';
        ax(3,jj).XAxisLocation = 'top';
        ax(3,jj).XLim = ax(1,jj).XLim;
        ax(3,jj).XTick = ax(1,jj).XTick;
        ax(3,jj).XTickLabels = velLabel;
      	ax(3,jj).XLabel.String = 'Peak Velocity (°/s)';
        ax(3,jj).XLabel.FontSize = ax(1,jj).YLabel.FontSize;
        
  	pp = pp + 1;
end
set(ax,'FontSize', 8, 'LineWidth', 1.5, 'XLim',[0.1 12])
set(ax, 'XTick', sort([ax(1,jj).XLim(1),2:2:12]))
linkaxes(ax,'x')
end