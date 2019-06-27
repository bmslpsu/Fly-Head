function [] = MakeFig_Sine_HeadFree_Grand_BODE()
%% MakeFig_Sine_HeadFree_Grand_BODE:
%   INPUTS:
%       -
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

[FILES,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','on');
FILES = cellstr(FILES)';

nAmp = length(FILES);
Amp = nan(nAmp,1);
for ww = 1:nAmp
    filedata = textscan(FILES{ww}, '%s', 'delimiter', '_');
    Amp(ww) = str2double(filedata{1}{3});
end

HeadFree = cell(nAmp,1);
for ww = 1:nAmp
    HeadFree{ww} = load(fullfile(root,FILES{ww}),'TRIAL','GRAND','U','N');
end

clearvars -except HeadFree Amp nAmp

%%
filename = 'MakeFig_Sine_HeadFree_Grand_BODE';

hold on
catIdx = 2;
xIdx = 1;
figNum = 1;

FIG = figure (figNum); clf
FIG.Color = 'w';
FIG.Position = [100 100 800 700];
movegui(FIG,'center')
FIG.Name = filename;
for ww = 1:nAmp
   FIG.Name = [FIG.Name '_' num2str(Amp(ww))];  
end

Freq = HeadFree{1}.U.freq{1};
nFreq = HeadFree{1}.N.freq;
GAIN = nan(nAmp,nFreq);
PHASE = nan(nAmp,nFreq);
pat2head = cell(6,1);

for ww = 1:nAmp
    for jj = 1:nFreq
        IOFreq = HeadFree{ww}.U.freq{1}(jj);
        
        pat.time    = HeadFree{ww}.TRIAL{2,jj}{1,1}.Time;
        pat.pos     = HeadFree{ww}.TRIAL{2,jj}{1,1}.X(:,xIdx);
        Pat = Fly(pat.pos,pat.time,[],IOFreq);
        
        head.time  	= HeadFree{ww}.GRAND{jj,catIdx}.Mean{1}{5};
        head.pos 	= HeadFree{ww}.GRAND{jj,catIdx}.Mean{1}{6}(:,xIdx);
        Head = Fly(head.pos,head.time,[],IOFreq);
        
        pat2head{jj} = IO_Class(Pat,Head);
        
        GAIN(ww,jj)     = pat2head{jj}.IOBodeGain(:,xIdx);
        PHASE(ww,jj)    = rad2deg(pat2head{jj}.IOBodePhaseDiff(:,xIdx));
        
        if PHASE(ww,jj)>10 && jj>=4
            PHASE(ww,jj) = PHASE(ww,jj) - 360;
        end
        
    end
    
    ax.gain = subplot(2,1,1); hold on
        ax.gain.FontSize = 12;
        ax.gain.YLabel.String = 'Gain';
        ax.gain.YLabel.FontSize = 14;
        ax.gain.XLabel.String = 'Frequency (Hz)';
        ax.gain.XLabel.FontSize = ax.gain.YLabel.FontSize;
        ax.gain.XLim = [0 12.5];
        ax.gain.YLim = [0 1];
        
        plot(Freq,GAIN(ww,:),'-o','LineWidth',2)
        
 	ax.phase = subplot(2,1,2); hold on
        ax.phase.FontSize = 12;
        ax.phase.YLabel.String = 'Phase';
        ax.phase.YLabel.FontSize = 14;
        ax.phase.XLabel.String = 'Frequency (Hz)';
        ax.phase.XLabel.FontSize = ax.phase.YLabel.FontSize;
        ax.phase.XLim = [0 12.5];
%         ax.phase.YLim = 180*[-1 1];
        
        plot(Freq,PHASE(ww,:),'-o','LineWidth',2)
end





end