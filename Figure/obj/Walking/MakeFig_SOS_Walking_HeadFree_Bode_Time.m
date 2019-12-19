function [FIG] = MakeFig_SOS_Walking_HeadFree_Bode_Time()
%% MakeFig_SOS_HeadFree_pat_head_Time:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'S:\Public\Audrey';

[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select Walking Data', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

Walking = load(fullfile(root,FILE{1}),'TRIAL','GRAND', 'FLY','U','N', 'T');

root = 'S:\Public\Audrey';
[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select Flight Data', root, 'MultiSelect','off');
FILE = cellstr(FILE)';
HeadFree = load(fullfile(root,FILE{1}),'TRIAL','GRAND', 'FLY','U','N', 'T');

%%
figNum = 1;
xIdx = 1;

cwalk = [0 0.4470 0.7410];
cfly = [0.8500 0.3250 0.0980];

filename = 'SOS_Walking_all_head_Time';

FIG = figure (figNum); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 3 1.5*3];
movegui(FIG,'center')
FIG.Name = filename;

ax.L = gca; hold on
ax.L.YColor = 'k';
ax.L.FontSize = 8;
ax.L.YLabel.String = ['Magnitude'];
ax.L.YLabel.FontSize = 10;
% ax.L.YLim =[0 10];
ax.L.YTick = [-270 0 180];
ax.L.XLabel.String = 'Frequency (Hz)';
ax.L.XLabel.Color = 'k';
ax.L.XLabel.FontSize = 10;
ax.L.XLim = [.5 10];
ax.L.XTick = [1, 3.1, 5.3, 7.4, 9.6];



plot(ax.L.XLim,[0 0],'--k','LineWidth',1);

% Grand Stats
ax1 = subplot(3,1,1) ; hold on
    ax1.FontSize = 8;
    ax1.XLim = [0 9.6]; 
    ax1.XTick = [1, 3.1, 5.3, 7.4, 9.6];
    ax1.YLim = [0 1];
	ax1.XLabel.String = 'Frequency (Hz)';
    ax1.XLabel.FontSize = 10;
 	ax1.YLabel.String = ['Gain (' char(176) '/' char(176) ')'];
    ax1.YLabel.FontSize = ax1.XLabel.FontSize;
    
    catIdx = 4;
    [~, WalkMag] = PlotPatch(Walking.GRAND{1,catIdx}.Mean{2}{5}(:,xIdx),Walking.GRAND{1,catIdx}.STD{2}{5}(:,xIdx),...
        Walking.GRAND{1,catIdx}.Mean{2}{4},1,Walking.N{1,1},cwalk,cwalk,0.5,1);
    
    catIdx = 5;
    [~, FlyMag] = PlotPatch(HeadFree.GRAND{1,catIdx}.Mean{2}{5}(:,xIdx),HeadFree.GRAND{1,catIdx}.STD{2}{5}(:,xIdx),...
        HeadFree.GRAND{1,catIdx}.Mean{2}{4},1, HeadFree.N{1,1},cfly,cfly,0.5,1);
    leg=legend([WalkMag FlyMag], 'Walking Fly', 'Flying Fly')
    leg.Box = 'off'
    %title('Bode Plot of Fly Response')
    
ax2 = subplot(3,1,2) ; hold on
    ax2.FontSize = ax1.FontSize;
    ax2.XLim = ax1.XLim;
    ax2.XTick = ax1.XTick;
    ax2.YLim = [-240 60];
   	ax2.XLabel.String = 'Frequency (Hz)';
    ax2.XLabel.FontSize = 10;
    ax2.XLabel.Color = 'k';
 	ax2.YLabel.String = ['Phase Difference(' char(176) ')'];
    ax2.YLabel.FontSize = ax1.XLabel.FontSize;
    ax2.YTick = -240:60:60;
    
    catIdx = 4;
    gphase = rad2deg(Walking.GRAND{1,catIdx}.CircMean{7}{6}(:,xIdx));
    %gphase(4) = gphase(4)-360;
    gphase(5) = gphase(5)-360;
    gstdphase = rad2deg(Walking.GRAND{1,catIdx}.CircSTD{7}{6}(:,xIdx));
    [~, WalkPhase] = PlotPatch(gphase,gstdphase,Walking.GRAND{1,catIdx}.Mean{2}{4},1,Walking.N{1,1},cwalk,cwalk,0.5,1);

    catIdx = 5;
    [~, FlyPhase] = PlotPatch(rad2deg(HeadFree.GRAND{1,catIdx}.CircMean{7}{6}(:,xIdx)),rad2deg(HeadFree.GRAND{1,catIdx}.CircSTD{7}{6}(:,xIdx)),...
        HeadFree.GRAND{1,catIdx}.Mean{2}{4},1,HeadFree.N{1,1},cfly,cfly,0.5,1);
    plot(ax1.XLim,[0 0],'--k','LineWidth',1);
 
ax3 = subplot(3,1,3) ; hold on
    ax3.FontSize = ax1.FontSize;
    ax3.XLim = ax1.XLim;
    ax3.XTick = ax1.XTick;
    %ax3.YLim = [-240 60];
   	ax3.XLabel.String = 'Frequency (Hz)';
    ax3.XLabel.FontSize = 10;
    ax3.XLabel.Color = 'k';
 	ax3.YLabel.String = ['Time Difference(ms)'];
    ax3.YLabel.FontSize = ax1.XLabel.FontSize;
    %ax3.YTick = -240:60:60;
    
    catIdx = 4;
    gphase = rad2deg(Walking.GRAND{1,catIdx}.CircMean{7}{6}(:,xIdx));
    %gphase(4) = gphase(4)-360;
    gphase(5) = gphase(5)-360;
    gtime = gphase/360.*((Walking.GRAND{1,catIdx}.Mean{2}{4}).^-1);
    gstdphase = rad2deg(Walking.GRAND{1,catIdx}.CircSTD{7}{6}(:,xIdx));
    gstdtime = gstdphase/360.*((Walking.GRAND{1,catIdx}.Mean{2}{4}).^-1)
    [~, WalkPhase] = PlotPatch(gtime*1000,gstdtime*1000,Walking.GRAND{1,catIdx}.Mean{2}{4},1,Walking.N{1,1},cwalk,cwalk,0.5,1);

    catIdx = 5;
    fphase = rad2deg(HeadFree.GRAND{1,catIdx}.CircMean{7}{6}(:,xIdx));
    ftime = fphase/360.*((HeadFree.GRAND{1,catIdx}.Mean{2}{4}).^-1);
    fstdphase = rad2deg(HeadFree.GRAND{1,catIdx}.CircSTD{7}{6}(:,xIdx));
    fstdtime = fstdphase/360.*((HeadFree.GRAND{1,catIdx}.Mean{2}{4}).^-1);
    [~, FlyPhase] = PlotPatch(ftime*1000,fstdtime*1000,...
        HeadFree.GRAND{1,catIdx}.Mean{2}{4},1,HeadFree.N{1,1},cfly,cfly,0.5,1);
    %plot(ax1.XLim,[0 0],'--k','LineWidth',1);
    

end