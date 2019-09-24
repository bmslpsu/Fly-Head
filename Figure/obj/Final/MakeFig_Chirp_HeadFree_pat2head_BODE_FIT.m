function [FIG] = MakeFig_Chirp_HeadFree_pat2head_BODE_FIT()
%% MakeFig_Chirp_HeadFree_pat2head_BODE_FIT: BODE head position for SOS
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';
% root = 'Q:\';
% Select files
[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

HeadFree = load(fullfile(root,FILE{1}),'TRIAL','GRAND','U','N');
%%
figNum = 1;
catIdx = 9;
xIdx = 1;

FIG = figure (figNum); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 4 4];
movegui(FIG,'center')
hold on

ax1 = subplot(2,1,1) ; hold on
    ax1.FontSize = 8;
    ax1.XLim = [0 10];
    ax1.YLim = [0 0.3];
    ax1.YLim = [0 1];
    ax1.XTickLabel = '';
    ax1.XLabel.FontSize = 8;
    ax1.XLabel.Color = 'none';
 	ax1.YLabel.String = ['Gain (' char(176) '/' char(176) ')'];
    ax1.YLabel.FontSize = ax1.XLabel.FontSize;
        
    [~,h.gain] = PlotPatch(HeadFree.GRAND{1,catIdx}.Mean{2}{5}(:,xIdx), HeadFree.GRAND{1,catIdx}.STD{2}{5}(:,xIdx), ...
        HeadFree.GRAND{1,catIdx}.Mean{2}{4}, 3, HeadFree.N{1,1}, 'b', [0.4 0.4 0.6], 0.5, 2);
	h.gain.Marker = '.';
    h.gain.MarkerSize = 20;
    
ax2 = subplot(2,1,2) ; hold on
    ax2.FontSize = ax1.FontSize;
    ax2.XLim = ax1.XLim;
    ax2.YLim = [-180 60];
   	ax2.XLabel.String = 'Frequency (Hz)';
    ax2.XLabel.FontSize = 8;
    ax2.XLabel.Color = 'k';
 	ax2.YLabel.String = ['Phase (' char(176) ')'];
    ax2.YLabel.FontSize = ax1.XLabel.FontSize;
    ax2.YTick = -300:60:60;
    
    [~,h.phase] = PlotPatch(rad2deg(HeadFree.GRAND{1,catIdx}.CircMean{7}{6}(:,xIdx)), rad2deg(HeadFree.GRAND{1,catIdx}.CircSTD{7}{6}(:,xIdx)), ...
        HeadFree.GRAND{1,catIdx}.Mean{2}{4}, 3, HeadFree.N{1,1}, 'b', [0.4 0.4 0.6], 0.5, 2);
    h.phase.Marker = '.';
    h.phase.MarkerSize = 20;
    
    plot(ax1.XLim,[0 0],'--k','LineWidth',1);
%%
% Freq = HeadFree.GRAND{1,catIdx}.Mean{2}{4};
% Gain = HeadFree.GRAND{1,catIdx}.Mean{2}{5}(:,xIdx);
% Phase = rad2deg(HeadFree.GRAND{1,catIdx}.CircMean{7}{6}(:,xIdx));

refIdx  = 1;
% headIdx = 2;
errIdx  = 4;
ampIdx  = 3;

Ts = HeadFree.GRAND{ampIdx,1}.Mean{2}{3};
% IN = HeadFree.GRAND{ampIdx,refIdx}.Mean{2}{6}(:,1);
U = HeadFree.TRIAL{1,ampIdx}{1,refIdx}.X(:,1);
Y = HeadFree.GRAND{ampIdx,errIdx}.Mean{2}{6}(:,1);
%%
mydata = iddata(Y,U,Ts);
                                    
% Transfer function estimation                              
 Options = tfestOptions;
 Options.Display = 'on';
 Options.WeightingFilter = [];
 np = 1;
 nz = 1;
 num = arrayfun(@(x)NaN(1,x), nz+1,'UniformOutput',false);
 den = arrayfun(@(x)[1, NaN(1,x)],np,'UniformOutput',false);

 % Prepare input/output delay                               
 iodValue = 0.01;
 iodFree = false;
 iodMin = 0.005;
 iodMax = 0.05;
 sysinit = idtf(num, den, 0);
 iod = sysinit.Structure.ioDelay;
 iod.Value = iodValue;
 iod.Free = iodFree;
 iod.Maximum = iodMax;
 iod.Minimum = iodMin;
 sysinit.Structure.ioDelay = iod;

 % Perform estimation using "sysinit" as template           
 chirp_tf = tfest(mydata, sysinit, Options);

end