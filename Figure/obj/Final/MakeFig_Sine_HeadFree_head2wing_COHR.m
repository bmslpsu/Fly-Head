function [FIG] = MakeFig_Sine_HeadFree_head2wing_COHR()
%% MakeFig_Sine_HeadFree_head2wing_COHR:
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
%%
figNum = 1;
filename = 'Sine_HeadFree_pat2head_COHR';
head2wing = 7;
xIdx = 1;
CC = prism(HeadFree{1}.N{1,3});

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
            FREQ.FlyMean{ww}(jj,kk)     = HeadFree{ww}.FLY{jj}{kk,head2wing}.Mean{4};
            COHR.FlyMean{ww}(jj,kk)     = HeadFree{ww}.FLY{jj}{kk,head2wing}.Mean{9}(:,xIdx);
            COHR.FlySTD{ww}(jj,kk)      = HeadFree{ww}.FLY{jj}{kk,head2wing}.STD{9}(:,xIdx);
        end
        FREQ.GrandMean{ww}(jj,1) 	= HeadFree{ww}.GRAND{jj,head2wing}.Mean{2}{4};
        COHR.GrandMean{ww}(jj,1)   	= HeadFree{ww}.GRAND{jj,head2wing}.Mean{2}{9}(:,xIdx);
        COHR.GrandSTD{ww}(jj,1)   	= HeadFree{ww}.GRAND{jj,head2wing}.STD{2}{9}(:,xIdx);
    end
end

FIG = figure (figNum) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 4 3];
FIG.Name = filename;
movegui(FIG,'center')

for ww = 1:nAmp
   FIG.Name = [FIG.Name '_' num2str(Amp(ww))];
end

hold on
ax = gca;
ax.FontSize = 8;
% ax.Title.String = [num2str(Amp) , char(176)];
ax.YLim = [0 1];
ax.YLabel.String = 'Coherence';
ax.YLabel.FontSize = 8;
ax.XLim = [0 13];
ax.XLabel.String = 'Frequency (Hz)';
ax.XLabel.FontSize = ax.YLabel.FontSize;
clear h
for ww = 1
    for jj = 1:HeadFree{ww}.N {1,3}
        [~,h(jj)] = PlotPatch(HeadFree{ww}.GRAND{jj,head2wing}.Mean{1}{7}(:,xIdx), HeadFree{ww}.GRAND{jj,head2wing}.STD{1}{7}(:,xIdx),...
            HeadFree{ww}.GRAND{jj,head2wing}.Mean{1}{8}, 2, HeadFree{ww}.N{1,1}, CC(jj,:), [0.4 0.4 0.6], 0.5, 2);
    end
end
uistack(h,'top')
end