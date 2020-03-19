function [FIG] = MakeFig_ChirpLog_HeadFree_head2wing_BODE_ALL_new()
%% MakeFig_ChirpLog_HeadFree_head2wing_BODE_ALL_new:
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

%%
clearvars -except HeadFree

filename = 'ChirpLog_HeadFree_head2wing_BODE_ALL_new';
catIdx = 7;
xIdx = 1;
CC = 'k';

% ALL Amp
FIG = figure (1); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 4*2.5 5];
FIG.Name = filename;
movegui(FIG,'center')
hold on

pp = 1;
for jj = 1:HeadFree.N{1,3}
	AMP(:,jj)   = HeadFree.U{1,3}{1}(jj);
    FREQ(:,jj)  = HeadFree.GRAND{jj,catIdx}.Mean{2}{1}(:,xIdx);
    VEL(:,jj) 	= AMP(:,jj)*2*pi*FREQ(:,jj);
    GAIN(:,jj)  = HeadFree.GRAND{jj,catIdx}.Mean{2}{2}(:,xIdx);
    PHASE(:,jj) = rad2deg(HeadFree.GRAND{jj,catIdx}.CircMean{7}{3}(:,xIdx));
    GSTD(:,jj)  = HeadFree.GRAND{jj,catIdx}.STD{2}{2}(:,xIdx);
    PSTD(:,jj)  = rad2deg(HeadFree.GRAND{jj,catIdx}.CircSTD{7}{3}(:,xIdx));
    
%     [b,a] = butter(2,0.5,'low');
%     [bb,aa] = butter(2,0.5,'low');
%     mff = 2;
%     PHASE(:,jj) = filtfilt(bb,aa,medfilt1(filtfilt(b,a,PHASE(:,jj)),mff));
%     GAIN(:,jj)  = filtfilt(bb,aa,medfilt1(filtfilt(b,a,GAIN(:,jj)),mff));
%     GSTD(:,jj)  = filtfilt(bb,aa,medfilt1(filtfilt(b,a,GSTD(:,jj)),mff));
%     PSTD(:,jj)  = filtfilt(bb,aa,medfilt1(filtfilt(b,a,PSTD(:,jj)),mff));
    
    ax1 = subplot(2,HeadFree.N{1,3},pp);
        hold on
        ax1.Title.String = [num2str(AMP(:,jj)) , char(176)];
        ax1.Title.Color = 'w';
        ax1.Title.FontSize = 16;
        ax1.FontSize = 12;
        ax1.YLabel.String = ['Gain (V/' char(176) ')'];
        ax1.YLabel.FontSize = 14;
        ax1.YLim = [0 0.7];
     	ax1.XLabel.String = 'Frequency (Hz)';
        ax1.XLabel.FontSize = ax1.YLabel.FontSize;
        ax1.XLabel.Color = 'w';
       	ax1.XLim = [0.1 12];
        ax1.XTickLabels = '';
        
        if pp>1
            ax1.YLabel.String = '';
            ax1.YTickLabels = '';
        end
        hold on

        h.patch = PlotPatch(GAIN(:,jj), GSTD(:,jj), FREQ(:,jj) ,...
                            2,HeadFree.N{1,1},CC,[0.4 0.4 0.6],0.5,1);
                        
        plot([0 12],[1 1],'--g','LineWidth',2);           
              
        ax1.XTick = unique(sort([min(ax1.XLim) ax1.XTick]));
        vel = round(AMP(:,jj)*2*pi*ax1.XTick);
        velLabel = cellfun(@(x) num2str(x), num2cell(vel), 'UniformOutput', false);
               
    ax2 = subplot(2,HeadFree.N{1,3},pp + HeadFree.N{1,3});
        hold on
        ax2.Title.String = [num2str(AMP(:,jj)) , char(176)];
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
        
        if pp>1
            ax2.YLabel.String = '';
            ax2.YTickLabels = '';
        end
        
        h.patch = PlotPatch(PHASE(:,jj), PSTD(:,jj), FREQ(:,jj) , 3,...
            HeadFree.N{1,1},CC,[0.4 0.4 0.6],0.5,1);
                
        plot([0 12],[0 0],'--g','LineWidth',2);
        
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
        
        clear ax
        pp = pp + 1;
end

%% One Amp
FIG = figure (1); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 4 4];
FIG.Name = filename;
movegui(FIG,'center')
hold on
CC = 'k';
amp = 3;

pp = 1;
for jj = amp
	AMP(:,jj)   = HeadFree.U{1,3}{1}(jj);
    FREQ(:,jj)  = HeadFree.GRAND{jj,catIdx}.Mean{2}{1}(:,xIdx);
    VEL(:,jj) 	= AMP(:,jj)*2*pi*FREQ(:,jj);
    GAIN(:,jj)  = HeadFree.GRAND{jj,catIdx}.Mean{2}{2}(:,xIdx);
    PHASE(:,jj) = rad2deg(HeadFree.GRAND{jj,catIdx}.CircMean{7}{3}(:,xIdx));
    GSTD(:,jj)  = HeadFree.GRAND{jj,catIdx}.STD{2}{2}(:,xIdx);
    PSTD(:,jj)  = rad2deg(HeadFree.GRAND{jj,catIdx}.CircSTD{7}{3}(:,xIdx));
    TD = 1000*(PHASE(:,jj)./360).*(1./FREQ(:,jj));
    TD_STD = 1000*(PSTD(:,jj)./360).*(1./FREQ(:,jj));
    TD(TD==inf) = 0;
    TD_STD(TD_STD==inf) = 0;
    
    [b,a] = butter(2,0.5,'low');
    [bb,aa] = butter(2,0.5,'low');
    mff = 2;
    PHASE(:,jj) = filtfilt(bb,aa,medfilt1(filtfilt(b,a,PHASE(:,jj)),mff));
    GAIN(:,jj)  = filtfilt(bb,aa,medfilt1(filtfilt(b,a,GAIN(:,jj)),mff));
    GSTD(:,jj)  = filtfilt(bb,aa,medfilt1(filtfilt(b,a,GSTD(:,jj)),mff));
    PSTD(:,jj)  = filtfilt(bb,aa,medfilt1(filtfilt(b,a,PSTD(:,jj)),mff));
    
    ax1 = subplot(3,1,pp);
        hold on
        ax1.Title.String = [num2str(AMP(:,jj)) , char(176)];
        ax1.Title.Color = 'w';
        ax1.Title.FontSize = 8;
        ax1.FontSize = 8;
        ax1.YLabel.String = ['Gain (V/' char(176) ')'];
        ax1.YLabel.FontSize = 8;
        ax1.YLim = [0 0.6];
     	ax1.XLabel.String = 'Frequency (Hz)';
        ax1.XLabel.FontSize = ax1.YLabel.FontSize;
        ax1.XLabel.Color = 'w';
       	ax1.XLim = [0.5 10];
        ax1.XTickLabels = '';
        
        if pp>1
            ax1.YLabel.String = '';
            ax1.YTickLabels = '';
        end
        hold on

        h.patch = PlotPatch(GAIN(:,jj), GSTD(:,jj), FREQ(:,jj) ,1 , HeadFree.N{1,1},CC,[0.4 0.4 0.6],0.5,2);
                        
        plot([0 12],[1 1],'--g','LineWidth',2);           
              
        ax1.XTick = unique(sort([min(ax1.XLim) ax1.XTick]));
        vel = round(AMP(:,jj)*2*pi*ax1.XTick);
        velLabel = cellfun(@(x) num2str(x), num2cell(vel), 'UniformOutput', false);
               
    ax2 = subplot(3,1,2);
        hold on
        ax2.Title.Color = 'k';
        ax2.Title.FontSize = ax1.Title.FontSize;
        ax2.FontSize = ax1.FontSize;
        ax2.YLabel.String = ['Phase Difference (' char(176) ')'];
        ax2.YLabel.FontSize = ax1.YLabel.FontSize;
        % ax2.XLabel.String = 'Frequency (Hz)';
        ax2.XLabel.FontSize = ax1.YLabel.FontSize;
      	ax2.XLabel.Color = 'k';
        ax2.YLim = rad2deg(pi*[-1 0.1]);
       	ax2.XLim = ax1.XLim;
        ax2.YTick = -180:60:180;
        
        if pp>1
            ax2.YLabel.String = '';
            ax2.YTickLabels = '';
        end
        
        h.patch = PlotPatch(PHASE(:,jj), PSTD(:,jj), FREQ(:,jj), 1, HeadFree.N{1,1},CC,[0.4 0.4 0.6],0.5,2);

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
        
    ax4 = subplot(3,1,3);
        hold on
        ax4.Title.Color = 'k';
        ax4.Title.FontSize = ax1.Title.FontSize;
        ax4.FontSize = ax1.FontSize;
        ax4.YLabel.String = 'Time Difference (ms)';
        ax4.YLabel.FontSize = ax1.YLabel.FontSize;
        ax4.XLabel.String = 'Frequency (Hz)';
        ax4.XLabel.FontSize = ax1.YLabel.FontSize;
      	ax4.XLabel.Color = 'k';
        ax4.YLim = [-180 15];
        ax4.YTick = -180:60:60;
       	ax4.XLim = ax1.XLim;
        ax4.XTick = ax3.XTick;
        
        h.patch = PlotPatch(TD, TD_STD, FREQ(:,jj), 1, HeadFree.N{1,1}, CC, [0.4 0.4 0.6], 0.5, 2);
                
        clear ax
        pp = pp + 1;
end

end