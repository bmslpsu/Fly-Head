function [FIG] = MakeFig_Chirp_HeadFreeadFixed_pat2head_BODE()
%% MakeFig_Chirp_HeadFreeadFixed_pat2head_BODE:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     : figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'F:\DATA\Rigid_Data\';

% Select head free chirp files
[FREE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head free chirp file', root, 'MultiSelect','on');
FREE = cellstr(FREE)';

% Select head fixed chirp files
[FIX,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head free chirp file', root, 'MultiSelect','on');
FIX = cellstr(FIX)';


HeadFree    = load(fullfile(root,FREE{1}),'TRIAL','FLY','GRAND','U','N');
HeadFixed   = load(fullfile(root,FIX{1}),'TRIAL','FLY','GRAND','U','N');

figNum = 1;
filename = 'Chirp_HeadFreeadFixed_pat2head_BODE'; % name of figure to save
catIdx = [8 3];
AmpIdx = 3;
xIdx = 1;

FIG = figure (figNum); clf % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 680 680];
FIG.Name = filename;
movegui(FIG,'center')
hold on


% Grand Stats
ax1 = subplot(2,1,1) ; hold on ; xlim([0.1 12]) ; ylim(0.3*[0 1]) ; title([num2str(HeadFree.U{1,3}{1}(AmpIdx)) , char(176)])
    ax1.FontSize = 12;
    ylabel(['Gain (V/' char(176) ')'],'FontSize',14)

    h.patch = PlotPatch(HeadFree.GRAND{AmpIdx,catIdx(1)}.Mean{2}{2}(:,xIdx),HeadFree.GRAND{AmpIdx,catIdx(1)}.STD{2}{2}(:,xIdx),...
HeadFree.GRAND{AmpIdx,catIdx(1)}.Mean{2}{1}(:,xIdx),3,HeadFree.N{1,1},'r',[0.4 0.4 0.6],0.5,2);
    
    h.patch = PlotPatch(HeadFixed.GRAND{AmpIdx,catIdx(2)}.Mean{2}{2}(:,xIdx),HeadFixed.GRAND{AmpIdx,catIdx(2)}.STD{2}{2}(:,xIdx),...
HeadFixed.GRAND{AmpIdx,catIdx(2)}.Mean{2}{1}(:,xIdx),3,HeadFixed.N{1,1},'c',[0.4 0.4 0.6],0.5,2);

    ax1.XTick = sort([0.1 ax1.XTick]);

ax2 = subplot(2,1,2) ; hold on ; xlim([0.1 12]) ; ylim(pi*[-1 1])
    h.Fly.Color(4) = 0.5;
    ax2.FontSize = 12;
    ylabel('Phase Difference (rad)','FontSize',14)
    xlabel('Frequency (Hz)','FontSize',14)

    h.patch = PlotPatch(HeadFree.GRAND{AmpIdx,catIdx(1)}.CircMean{9}{3}(:,xIdx),HeadFree.GRAND{AmpIdx,catIdx(1)}.CircSTD{9}{3}(:,xIdx),...
HeadFree.GRAND{AmpIdx,catIdx(1)}.Mean{2}{1}(:,xIdx),3,HeadFree.N{1,1},'r',[0.4 0.4 0.6],0.5,2);

    h.patch = PlotPatch(medfilt1(HeadFixed.GRAND{AmpIdx,catIdx(2)}.CircMean{9}{3}(:,xIdx),1),HeadFixed.GRAND{AmpIdx,catIdx(2)}.CircSTD{9}{3}(:,xIdx),...
HeadFixed.GRAND{AmpIdx,catIdx(2)}.Mean{2}{1}(:,xIdx),3,HeadFixed.N{1,1},'c',[0.4 0.4 0.6],0.5,2);

    plot([0 12],[0 0],'--g','LineWidth',2);

    ax2.XTick = ax1.XTick;
        
end