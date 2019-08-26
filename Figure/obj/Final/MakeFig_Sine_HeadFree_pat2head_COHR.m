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
%%
figNum = 1;
filename = 'Sine_HeadFree_COHR';
headIdx = 5;
wingIdx = 8;
xIdx = 1;

% Fly Stats
FREQ.GrandMean	= cell(nAmp,2);
COHR.GrandMean 	= cell(nAmp,2);
COHR.GrandSTD	= cell(nAmp,2);
for ww = 1:nAmp
    for jj = 1:HeadFree{ww}.N{1,3}
        FREQ.GrandMean{ww,1}(jj,1)      = HeadFree{ww}.GRAND{jj,headIdx}.Mean{2}{4};
        COHR.GrandMean{ww,1}(jj,1)   	= HeadFree{ww}.GRAND{jj,headIdx}.Mean{2}{9}(:,xIdx);
        COHR.GrandSTD{ww,1}(jj,1)   	= HeadFree{ww}.GRAND{jj,headIdx}.STD{2}{9}(:,xIdx);
        
     	FREQ.GrandMean{ww,2}(jj,1)      = HeadFree{ww}.GRAND{jj,wingIdx}.Mean{2}{4};
        COHR.GrandMean{ww,2}(jj,1)   	= HeadFree{ww}.GRAND{jj,wingIdx}.Mean{2}{9}(:,xIdx);
        COHR.GrandSTD{ww,2}(jj,1)   	= HeadFree{ww}.GRAND{jj,wingIdx}.STD{2}{9}(:,xIdx);
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
for ww = 1
    for jj = 1:HeadFree{ww}.N {1,3}
%         PlotPatch(HeadFree{ww}.GRAND{jj,headIdx}.Mean{1}{7}(:,xIdx), HeadFree{ww}.GRAND{jj,headIdx}.STD{1}{7}(:,xIdx),...
%             HeadFree{ww}.GRAND{jj,headIdx}.Mean{1}{8}, 2, HeadFree{ww}.N{1,1}, 'b', [0.4 0.4 0.6], 0.5, 1);
%         
%         PlotPatch(HeadFree{ww}.GRAND{jj,wingIdx}.Mean{1}{7}(:,xIdx), HeadFree{ww}.GRAND{jj,wingIdx}.STD{1}{7}(:,xIdx),...
%             HeadFree{ww}.GRAND{jj,wingIdx}.Mean{1}{8}, 2, HeadFree{ww}.N{1,1}, 'r', [0.4 0.4 0.6], 0.5, 1);
        
%         plot(HeadFree{ww}.GRAND{jj,catIdx}.Mean{1}{8},HeadFree{ww}.GRAND{jj,catIdx}.Mean{1}{7}(:,xIdx),'Color',CC(jj,:),...
%             'LineWidth',1)
    end
	errorbar(FREQ.GrandMean{ww,1},COHR.GrandMean{ww,1},2*COHR.GrandSTD{ww,1}./sqrt(HeadFree{ww}.N{1,1}),'-b','LineWidth',2)
 	errorbar(FREQ.GrandMean{ww,2},COHR.GrandMean{ww,2},2*COHR.GrandSTD{ww,2}./sqrt(HeadFree{ww}.N{1,1}),'-r','LineWidth',2)
end

end