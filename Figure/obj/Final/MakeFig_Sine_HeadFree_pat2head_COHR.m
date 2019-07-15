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

CC = prism(HeadFree{ww}.N {1,3});
for ww = 3
    for jj = 1:HeadFree{ww}.N {1,3}
        PlotPatch(HeadFree{ww}.GRAND{jj,5}.Mean{1}{7}(:,xIdx), HeadFree{ww}.GRAND{jj,5}.STD{1}{7}(:,xIdx),...
            HeadFree{ww}.GRAND{jj,5}.Mean{1}{8}, 2, HeadFree{ww}.N{1,1}, 'b', [0.4 0.4 0.6], 0.5, 1);
        
        PlotPatch(HeadFree{ww}.GRAND{jj,8}.Mean{1}{7}(:,xIdx), HeadFree{ww}.GRAND{jj,8}.STD{1}{7}(:,xIdx),...
            HeadFree{ww}.GRAND{jj,8}.Mean{1}{8}, 2, HeadFree{ww}.N{1,1}, 'r', [0.4 0.4 0.6], 0.5, 1);
%         plot(HeadFree{ww}.GRAND{jj,catIdx}.Mean{1}{8},HeadFree{ww}.GRAND{jj,catIdx}.Mean{1}{7}(:,xIdx),'Color',CC(jj,:),...
%             'LineWidth',1)
    end
% 	errorbar(FREQ.GrandMean{ww},COHR.GrandMean{ww},1*COHR.GrandSTD{ww},'-or','LineWidth',3)
end

end