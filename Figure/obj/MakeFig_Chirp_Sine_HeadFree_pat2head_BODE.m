function [FIG] = MakeFig_Chirp_Sine_HeadFree_pat2head_BODE()
%% MakeFig_Chirp_Sine_HeadFree_pat2head_BODE:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     : figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';
root = 'C:\Users\boc5244\Downloads\';

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

figNum = 4;
filename = 'Chirp_Sine_HeadFree_pat2head_BODE'; % name of figure to save
catIdx = 5;
xIdx = 1;
chirpIdx = 3;

FREQ.GrandMean	= cell(nAmp,1);
GAIN.GrandMean 	= cell(nAmp,1);
PHASE.GrandMean	= cell(nAmp,1);
GAIN.GrandSTD	= cell(nAmp,1);
PHASE.GrandSTD	= cell(nAmp,1);
for ww = 1:nAmp % amplitudes
    for jj = 1:HeadFree{ww}.N{1,3}-1 % frequencies
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
    ax1 = subplot(2,1,1);
        hold on
%         ax1.Title.String = [num2str(Amp) , char(176)];
        ax1.Title.FontSize = 16;
        ax1.FontSize = 12;
        ax1.YLabel.String = ['Gain (' char(176) '/' char(176) ')'];
        ax1.YLabel.FontSize = 14;
        ax1.YLim = [0 1];
       	ax1.XLim = [0.1 8];
        ax1.XTickLabels = '';
        
        h.patch = PlotPatch(HeadFree{nAmp+1}.GRAND{chirpIdx,catIdx}.Mean{2}{2}(:,xIdx),...
                            HeadFree{nAmp+1}.GRAND{chirpIdx,catIdx}.STD{2}{2}(:,xIdx),...
                            HeadFree{nAmp+1}.GRAND{chirpIdx,catIdx}.Mean{2}{1}(:,xIdx),...
                            3,HeadFree{nAmp+1}.N{1,1},'b',[0.4 0.4 0.6],0.5,2);
        
        h.Fly = errorbar(FREQ.GrandMean{ww},GAIN.GrandMean{ww},2*GAIN.GrandSTD{ww},'-ok','LineWidth',3);
        
        ax1.XTick = sort([0.1 ax1.XTick]);
        vel = round(Amp*2*pi*ax1.XTick);
        velLabel = cellfun(@(x) num2str(x), num2cell(vel), 'UniformOutput', false);
        
    ax2 = subplot(2,1,2) ; hold on;
        hold on
        ax2.Position = ax2.Position + [0 0.05 0 0];
%         ax2.Title.String = [num2str(Amp) , char(176)];
        ax2.Title.FontSize = 16;
        ax2.FontSize = ax1.FontSize;
        ax2.YLabel.String = ['Phase Difference (' char(176) ')'];
        ax2.YLabel.FontSize = ax1.YLabel.FontSize;
        ax2.XLabel.String = 'Frequency (Hz)';
        ax2.XLabel.FontSize = ax1.YLabel.FontSize;
        ax2.YLim = rad2deg(pi*[-1 1]);
       	ax2.XLim = ax1.XLim;
        
        h.patch = PlotPatch(rad2deg(HeadFree{nAmp+1}.GRAND{3,catIdx}.CircMean{9}{3}(:,xIdx)),...
            rad2deg(HeadFree{nAmp+1}.GRAND{3,catIdx}.CircSTD{9}{3}(:,xIdx)),...
            HeadFree{nAmp+1}.GRAND{3,catIdx}.Mean{2}{1}(:,xIdx),3,HeadFree{nAmp+1}.N{1,1},'b',[0.4 0.4 0.6],0.5,2);
        
        h.Fly = errorbar(FREQ.GrandMean{ww},rad2deg(PHASE.GrandMean{ww}),rad2deg(2*PHASE.GrandSTD{ww}),'-ok','LineWidth',3);
        
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
end

end