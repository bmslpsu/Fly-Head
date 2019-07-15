function [FIG] = MakeFig_ChirpLog_HeadFree_pat2head_COHR_all()
%% MakeFig_ChirpLog_HeadFree_pat2head_COHR_all:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     : figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

[CHIRP,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select chirp file', root, 'MultiSelect','on');
CHIRP = cellstr(CHIRP)';

HeadFree = load(fullfile(root,CHIRP{1}),'GRAND','U','N');
%%
figNum = 1;
filename = 'ChirpLog_HeadFree_pat2head_COHR_all';
catIdx = 5;
xIdx = 1;
CC = [0 0 0];
offset = 0.2;

FIG = figure (figNum); clf
FIG.Color = 'w';
FIG.Position = [100 100 1200 400];
FIG.Name = filename;
movegui(FIG,'center')
hold on

pp = 1;
for jj = 1:HeadFree.N{1,3}     
    Amp = HeadFree.U{1,3}{1}(jj);
    ax1 = subplot(1,HeadFree.N{1,3},pp);
        hold on
        initPos = ax1.Position;
        ax1.Position = ax1.Position + [0 offset 0 -offset];
%         ax1.Title.String = [num2str(Amp) , char(176)];
        ax1.Title.Color = 'k';
        ax1.Title.FontSize = 16;
        ax1.FontSize = 12;
        ax1.YLabel.String = 'Coherence';
        ax1.YLabel.FontSize = 14;
        ax1.YLim = [0 1.0];
        ax1.YTick = unique(sort([ax1.YTick ax1.YLim(2)]));
     	ax1.XLabel.String = 'Frequency (Hz)';
        ax1.XLabel.Color = 'k';
        ax1.XLabel.FontSize = ax1.YLabel.FontSize;
       	ax1.XLim = [0.3 10];
        
        if pp>1
            ax1.YLabel.String = '';
            ax1.YTickLabels = '';
        end
        hold on
    
        FREQ = HeadFree.GRAND{jj,catIdx}.Mean{1}{8};
     	COHR = HeadFree.GRAND{jj,catIdx}.Mean{1}{7}(:,xIdx);
     	STD  = HeadFree.GRAND{jj,catIdx}.STD{1}{7}(:,xIdx);
        
       	PlotPatch(COHR, STD, FREQ , 3, HeadFree.N{1,1}, CC, [0.4 0.4 0.6], 0.5, 2);
                                     
        ax1.XTick = unique(sort([min(ax1.XLim) ax1.XTick]));
        vel = round(Amp*2*pi*ax1.XTick);
        velLabel = cellfun(@(x) num2str(x), num2cell(vel), 'UniformOutput', false);               
        
       	ax3 = axes;
        ax3.Position = initPos + [0 0.05 0 0];
        ax3.FontSize = ax1.FontSize ;
        ax3.Color = 'none';
        ax3.YAxisLocation = 'right';
        ax3.YAxis.Color = 'none';
        ax3.XAxisLocation = 'bottom';
        ax3.XLim = ax1.XLim;
        ax3.XTick = ax1.XTick;
        ax3.XTickLabels = velLabel;
      	ax3.XLabel.String = ['Peak Velocity (' char(176) '/s)'];
        ax3.XLabel.FontSize = ax1.YLabel.FontSize;
        
        clear ax
        pp = pp + 1;
end

end