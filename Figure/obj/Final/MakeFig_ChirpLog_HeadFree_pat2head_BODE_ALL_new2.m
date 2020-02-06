function [FIG] = MakeFig_ChirpLog_HeadFree_pat2head_BODE_ALL_new2()
%% MakeFig_ChirpLog_HeadFree_pat2head_BODE_ALL_new:
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
filename = 'ChirpLog_HeadFree_pat2head_BODE_ALL_new';
catIdx = 6;
xIdx = 1;
CC = [0 0 0.7];
CC = [0.7 0 0];
y_lim = 0.2;

FIG = figure (1); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 10 5];
FIG.Name = filename;
movegui(FIG,'center')
hold on
Fly
AMP     = [];
FREQ    = [];
VEL     = [];
GAIN    = [];
PHASE   = [];
clear ax
pp = 1;
for jj = 1:HeadFree.N{1,3} % amplitudes
	AMP(:,jj)      = HeadFree.U{1,3}{1}(jj);
    FREQ(:,jj)     = HeadFree.GRAND{jj,catIdx}.Mean{2}{1}(:,xIdx);
    VEL(:,jj)      = AMP(:,jj)*2*pi*FREQ(:,jj);
    GAIN(:,jj)     = abs(HeadFree.GRAND{jj,catIdx}.Mean{1}{17}(:,xIdx));
    PHASE(:,jj)    = rad2deg(angle(HeadFree.GRAND{jj,catIdx}.Mean{1}{17}(:,xIdx)));
    GSTD(:,jj)     = abs(HeadFree.GRAND{jj,catIdx}.STD{1}{17}(:,xIdx));
    PSTD(:,jj)     = rad2deg((HeadFree.GRAND{jj,catIdx}.STD{1}{17}(:,xIdx)));
    
 	[b,a] = butter(2,0.5,'low');
    [bb,aa] = butter(2,0.5,'low');
    mff = 1;
    PHASE(:,jj) = filtfilt(bb,aa,medfilt1(filtfilt(b,a,PHASE(:,jj)),mff));
    GAIN(:,jj)  = filtfilt(bb,aa,medfilt1(filtfilt(b,a,GAIN(:,jj)),mff));
    GSTD(:,jj)  = filtfilt(bb,aa,medfilt1(filtfilt(b,a,GSTD(:,jj)),mff));
    PSTD(:,jj)  = filtfilt(bb,aa,medfilt1(filtfilt(b,a,PSTD(:,jj)),mff));
    
    ax(1) = subplot(2,HeadFree.N{1,3},pp);
        hold on
        ax(1).YLabel.String = ['Gain (' char(176) '/' char(176) ')'];
        ax(1).YLabel.FontSize = 10;
        ax(1).YLim = [0 y_lim];
        ax(1).YTick = unique(sort([ax(1).YTick ax(1).YLim(2)]));

        h.patch = PlotPatch(GAIN(:,jj), GSTD(:,jj), FREQ(:,jj) ,...
                            1,HeadFree.N{1,1},CC,[0.4 0.4 0.6],0.5,1);
                               
        ax(1).XTick = unique(sort([min(ax(1).XLim) 2:2:12]));
        vel = round(AMP(:,jj)*2*pi*ax(1).XTick);
        velLabel = cellfun(@(x) num2str(x), num2cell(vel), 'UniformOutput', false);
               
    ax(2) = subplot(2,HeadFree.N{1,3},pp + HeadFree.N{1,3});
        hold on
        ax(2).Title.String = [num2str(AMP(:,jj)) , char(176)];
        ax(2).Title.Color = 'k';
        ax(2).Title.FontSize = ax(1).Title.FontSize;
        ax(2).YLabel.String = ['Phase Difference (' char(176) ')'];
        ax(2).YLabel.FontSize = ax(1).YLabel.FontSize;
        ax(2).XLabel.String = 'Frequency (Hz)';
        ax(2).XLabel.FontSize = ax(1).YLabel.FontSize;
      	ax(2).XLabel.Color = 'k';
        ax(2).YLim = rad2deg(pi*[-1 1]);
        ax(2).YTick = -180:60:180;
        
        h.patch = PlotPatch(PHASE(:,jj), PSTD(:,jj),FREQ(:,jj),...
            2,HeadFree.N{1,1},CC,[0.4 0.4 0.6],0.5,1);
                
        plot([0 12],[0 0],'--k','LineWidth',1)
                
	ax(3) = axes;
        ax(3).Position = ax(1).Position + [0 -0.00 0 0];
        ax(3).FontSize = ax(1).FontSize ;
        ax(3).Color = 'none';
        ax(3).YAxisLocation = 'right';
        ax(3).YAxis.Color = 'none';
        ax(3).XAxisLocation = 'top';
        ax(3).XLim = ax(1).XLim;
        ax(3).XTick = ax(1).XTick;
        ax(3).XTickLabels = velLabel;
      	ax(3).XLabel.String = ['Peak Velocity (' char(176) '/s)'];
        ax(3).XLabel.FontSize = ax(1).YLabel.FontSize;
        
        pp = pp + 1;
      	set(ax,'FontSize', 8, 'XLim',[0.1 12])
        set(ax, 'XTick', sort([ax(1).XLim(1),2:2:12]))
        linkaxes(ax,'x')
end

end