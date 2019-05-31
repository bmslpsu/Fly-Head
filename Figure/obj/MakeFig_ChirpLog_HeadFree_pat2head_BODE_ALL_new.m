function [FIG] = MakeFig_ChirpLog_HeadFree_pat2head_BODE_ALL_new()
%% MakeFig_ChirpLog_HeadFree_pat2head_BODE_ALL_new:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     : figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

% Select chirp files
[CHIRP,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select chirp file', root, 'MultiSelect','on');
CHIRP = cellstr(CHIRP)';

HeadFree = load(fullfile(root,CHIRP{1}),'GRAND','U','N');

figNum = 1;
filename = 'ChirpLog_HeadFree_pat2head_BODE_ALL_new'; % name of figure to save
catIdx = 5;
xIdx = 1;
CC = [0 0 0.7];

FIG = figure (figNum); clf % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1200 650];
% FIG.OuterPosition = FIG.Position + [0 0 1 400];
FIG.Name = filename;
movegui(FIG,'center')
hold on

AMP     = [];
FREQ    = [];
VEL     = [];
GAIN    = [];
PHASE   = [];

pp = 1;
for jj = 1:HeadFree.N{1,3} % amplitudes
	AMP(:,jj)   = HeadFree.U{1,3}{1}(jj);
    FREQ(:,jj)  = HeadFree.GRAND{jj,catIdx}.Mean{2}{1}(:,xIdx);
    VEL(:,jj) 	= AMP(:,jj)*2*pi*FREQ(:,jj);
    GAIN(:,jj)  = HeadFree.GRAND{jj,catIdx}.Mean{2}{2}(:,xIdx);
    PHASE(:,jj) = rad2deg(HeadFree.GRAND{jj,catIdx}.CircMean{9}{3}(:,xIdx));
    
    ax1 = subplot(2,HeadFree.N{1,3},pp);
        hold on
        ax1.Title.String = [num2str(AMP(:,jj)) , char(176)];
        ax1.Title.Color = 'w';
        ax1.Title.FontSize = 16;
        ax1.FontSize = 12;
        ax1.YLabel.String = ['Gain (' char(176) '/' char(176) ')'];
        ax1.YLabel.FontSize = 14;
        ax1.YLim = [0 1.2];
     	ax1.XLabel.String = 'Frequency (Hz)';
        ax1.XLabel.FontSize = ax1.YLabel.FontSize;
        ax1.XLabel.Color = 'w';
       	ax1.XLim = [0.3 10];
        ax1.XTickLabels = '';
%         ax1.OuterPosition = ax1.OuterPosition + [0 0 0 0.1];
        
        if pp>1
            ax1.YLabel.String = '';
            ax1.YTickLabels = '';
        end
        hold on

        h.patch = PlotPatch(GAIN(:,jj),...
                            HeadFree.GRAND{jj,catIdx}.STD{2}{2}(:,xIdx),...
                            FREQ(:,jj) ,...
                            3,HeadFree.N{1,1},CC,[0.4 0.4 0.6],0.5,2);
                        
        plot([0 12],[1 1],'--g','LineWidth',2);           
              
        ax1.XTick = unique(sort([min(ax1.XLim) ax1.XTick]));
        vel = round(AMP(:,jj)*2*pi*ax1.XTick);
        velLabel = cellfun(@(x) num2str(x), num2cell(vel), 'UniformOutput', false);
               
    ax2 = subplot(2,HeadFree.N{1,3},pp + HeadFree.N{1,3});
        hold on
%         ax2.Position = ax2.Position + [0 0.05 0 0];
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
        
        if pp>1
            ax2.YLabel.String = '';
            ax2.YTickLabels = '';
        end
        
        h.patch = PlotPatch(PHASE(:,jj),...
            rad2deg(HeadFree.GRAND{jj,catIdx}.CircSTD{9}{3}(:,xIdx)),...
            FREQ(:,jj) ,3,HeadFree.N{1,1},CC,[0.4 0.4 0.6],0.5,2);
                
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
%%
ALL = sortrows([FREQ(:) , VEL(:), GAIN(:)],1);
ALL = ALL(1:962,:);

FIG2 = figure (2) ; clf
FIG2.Color = 'w';
hold on
CC = jet(4);
ax = gca;
for jj = 1:length(AMP)
%     hh = plot3(FREQ(:,jj),VEL(:,jj),GAIN(:,jj),'-','LineWidth',2) ; hold on
    
    hh =    PlotPatch(GAIN(:,jj),...
            HeadFree.GRAND{jj,catIdx}.STD{2}{2}(:,xIdx),...
            FREQ(:,jj) ,...
            2,HeadFree.N{1,1},CC(jj,:),[0.4 0.4 0.6],0.5,2);
    
%     hh.Color(4) = 0.7;
end
xlim([25 1200])
zlim([0 1.2])
ylim([0 1.2])
ax.XLabel.String = 'Frequency (Hz)';
ax.YLabel.String = 'Velocity';
ax.ZLabel.String = 'Gain';
leglabel = cellfun(@(x) num2str(x), num2cell(AMP),'UniformOutput',false);
% legend(leglabel)
% view(90,0)
% view(0,0)


end