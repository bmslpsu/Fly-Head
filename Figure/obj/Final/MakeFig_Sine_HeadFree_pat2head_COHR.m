function [FIG] = MakeFig_Sine_HeadFree_pat2head_COHR()
%% MakeFig_Sine_HeadFree_pat2head_COHR:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     : figure handle
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
    HeadFree{ww} = load(fullfile(root,FILES{ww}),'TRIAL','FLY','GRAND','U','N');
end

figNum = 1;
filename = 'Sine_HeadFree_pat2head_COHR';
catIdx = 5;
xIdx = 1;

% Fly Stats
FREQ.FlyMean    = cell(nAmp,1);
COHR.FlyMean  	= cell(nAmp,1);
FREQ.GrandMean	= cell(nAmp,1);
COHR.GrandMean 	= cell(nAmp,1);
COHR.FlySTD     = cell(nAmp,1);
COHR.GrandSTD	= cell(nAmp,1);
for ww = 1:nAmp
    for jj = 1:HeadFree{ww}.N{1,3}
        for kk = 1:HeadFree{ww}.N{1,1}
            FREQ.FlyMean{ww}(jj,kk)     = HeadFree{ww}.FLY{jj}{kk,catIdx}.Mean{4};
            COHR.FlyMean{ww}(jj,kk)     = HeadFree{ww}.FLY{jj}{kk,catIdx}.Mean{9}(:,xIdx);
            COHR.FlySTD{ww}(jj,kk)      = HeadFree{ww}.FLY{jj}{kk,catIdx}.STD{9}(:,xIdx);
        end
        FREQ.GrandMean{ww}(jj,1) 	= HeadFree{ww}.GRAND{jj,catIdx}.Mean{2}{4};
        COHR.GrandMean{ww}(jj,1)   	= HeadFree{ww}.GRAND{jj,catIdx}.Mean{2}{9}(:,xIdx);
        COHR.GrandSTD{ww}(jj,1)   	= HeadFree{ww}.GRAND{jj,catIdx}.STD{2}{9}(:,xIdx);
    end
end

FIG = figure (figNum) ; clf
FIG.Color = 'w';
FIG.Position = [100 100 680 500];
FIG.Name = filename;
movegui(FIG,'center')

for ww = 1:nAmp
   FIG.Name = [FIG.Name '_' num2str(Amp(ww))];  
end

hold on
ax = gca;
ax.FontSize = 12;
ax.Title.String = [num2str(Amp) , char(176)];
ax.YLim = [0 1];
ax.YLabel.String = 'Coherence';
ax.YLabel.FontSize = 14;
ax.XLim = [0 13];
ax.XLabel.String = 'Frequency (Hz)';
ax.XLabel.FontSize = ax.YLabel.FontSize;

CC = jet(HeadFree{ww}.N {1,3});
for ww = 1:nAmp
    for jj = 1:HeadFree{ww}.N {1,3}
        PlotPatch(HeadFree{ww}.GRAND{jj,catIdx}.Mean{1}{7}(:,xIdx), HeadFree{ww}.GRAND{jj,catIdx}.STD{1}{7}(:,xIdx),...
            HeadFree{ww}.GRAND{jj,catIdx}.Mean{1}{8}, 3, HeadFree{ww}.N{1,1}, CC(jj,:), [0.4 0.4 0.6], 0.5, 2);
    end
	errorbar(FREQ.GrandMean{ww},COHR.GrandMean{ww},1*COHR.GrandSTD{ww},'-or','LineWidth',3)
end

end