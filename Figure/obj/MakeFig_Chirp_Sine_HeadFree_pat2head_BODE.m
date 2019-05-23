function [FIG] = MakeFig_Chirp_Sine_HeadFree_pat2head_BODE()
%% MakeFig_Chirp_Sine_HeadFree_pat2head_BODE:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     : figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'F:\DATA\Rigid_Data\';

% Select sine files
[SINE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select sine files', root, 'MultiSelect','on');
SINE = cellstr(SINE)';

% Select chirp files
[CHIRP,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select chirp file', root, 'MultiSelect','on');
CHIRP = cellstr(CHIRP)';

nAmp = length(SINE);
Amp = nan(nAmp,1);
for ww = 1:nAmp
    filedata = textscan(SINE{ww}, '%s', 'delimiter', '_');
    Amp(ww) = str2double(filedata{1}{3});
end

HeadFree = cell(nAmp+1,1);
for ww = 1:nAmp
    HeadFree{ww} = load(fullfile(root,SINE{ww}),'TRIAL','FLY','GRAND','U','N');
end

HeadFree{nAmp+1} = load(fullfile(root,CHIRP{ww}),'TRIAL','FLY','GRAND','U','N');


figNum = 1;
filename = 'Chirp_Sine_HeadFree_pat2head_BODE'; % name of figure to save
catIdx = 5;
xIdx = 1;

FREQ.GrandMean	= cell(nAmp,1);
GAIN.GrandMean 	= cell(nAmp,1);
PHASE.GrandMean	= cell(nAmp,1);
GAIN.GrandSTD	= cell(nAmp,1);
PHASE.GrandSTD	= cell(nAmp,1);
for ww = 1:nAmp % amplitudes
    for jj = 1:HeadFree{ww}.N{1,3} % frequencies
        FREQ.GrandMean{ww}(jj,1) 	= HeadFree{ww}.GRAND{jj,catIdx}.Mean{2}{4};
        GAIN.GrandMean{ww}(jj,1)   	= HeadFree{ww}.GRAND{jj,catIdx}.Mean{2}{5}(:,xIdx);
        PHASE.GrandMean{ww}(jj,1) 	= HeadFree{ww}.GRAND{jj,catIdx}.CircMean{9}{6}(:,xIdx);
        GAIN.GrandSTD{ww}(jj,1)   	= HeadFree{ww}.GRAND{jj,catIdx}.STD{2}{5}(:,xIdx);
        PHASE.GrandSTD{ww}(jj,1) 	= HeadFree{ww}.GRAND{jj,catIdx}.CircSTD{9}{6}(:,xIdx);
    end
end

FIG = figure (figNum); clf % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 680 680];
FIG.Name = filename;
movegui(FIG,'center')
for ww = 1:nAmp
   FIG.Name = [FIG.Name '_' num2str(Amp(ww))];  
end
hold on


% Grand Stats
for ww = 1:nAmp % amplitudes
    ax1 = subplot(2,1,1) ; hold on ; xlim([0.1 12]) ; ylim(1*[0 1]) ; title([num2str(Amp) , char(176)])
        ax1.FontSize = 12;
        ylabel(['Gain (V/' char(176) ')'],'FontSize',14)
                
        h.patch = PlotPatch(HeadFree{nAmp+1}.GRAND{3,catIdx}.Mean{2}{2}(:,xIdx),HeadFree{nAmp+1}.GRAND{3,catIdx}.STD{2}{2}(:,xIdx),...
    HeadFree{nAmp+1}.GRAND{3,catIdx}.Mean{2}{1}(:,xIdx),3,HeadFree{nAmp+1}.N{1,1},'b',[0.4 0.4 0.6],0.5,2);
        
        h.Fly = errorbar(FREQ.GrandMean{ww},GAIN.GrandMean{ww},2*GAIN.GrandSTD{ww},'-ok','LineWidth',3);
        
        ax1.XTick = sort([0.1 ax1.XTick]);

    ax2 = subplot(2,1,2) ; hold on ; xlim([0.1 12]) ; ylim(pi*[-1 1])
        h.Fly.Color(4) = 0.5;
        ax2.FontSize = 12;
        ylabel('Phase Difference (rad)','FontSize',14)
        xlabel('Frequency (Hz)','FontSize',14)
        
        h.patch = PlotPatch(HeadFree{nAmp+1}.GRAND{3,catIdx}.CircMean{9}{3}(:,xIdx),HeadFree{nAmp+1}.GRAND{3,catIdx}.CircSTD{9}{3}(:,xIdx),...
    HeadFree{nAmp+1}.GRAND{3,catIdx}.Mean{2}{1}(:,xIdx),3,HeadFree{nAmp+1}.N{1,1},'b',[0.4 0.4 0.6],0.5,2);
        
        h.Fly = errorbar(FREQ.GrandMean{ww},PHASE.GrandMean{ww},2*PHASE.GrandSTD{ww},'-ok','LineWidth',3);
        
        plot([0 12],[0 0],'--g','LineWidth',2);
        
        ax2.XTick = ax1.XTick;
        
end

end