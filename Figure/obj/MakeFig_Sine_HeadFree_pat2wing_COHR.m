function [FIG] = MakeFig_Sine_HeadFree_pat2wing_COHR()
%% MakeFig_Sine_HeadFree_pat2wing_COHR:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     : figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';
figNum = 1;
filename = 'Sine_HeadFree_pat2wing_COHR'; % name of figure to save
catIdx = 6; % pat2head
xIdx = 1;

% Select files
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

% Fly Stats
FREQ.FlyMean    = cell(nAmp,1);
COHR.FlyMean  	= cell(nAmp,1);
FREQ.GrandMean	= cell(nAmp,1);
COHR.GrandMean 	= cell(nAmp,1);
COHR.FlySTD     = cell(nAmp,1);
COHR.GrandSTD	= cell(nAmp,1);
for ww = 1:nAmp % amplitudes
    for jj = 1:HeadFree{ww}.N{1,3} % frequencies
        for kk = 1:HeadFree{ww}.N{1,1} % flys
            FREQ.FlyMean{ww}(jj,kk)     = HeadFree{ww}.FLY{jj}{kk,catIdx}.Mean{4};
            COHR.FlyMean{ww}(jj,kk)     = HeadFree{ww}.FLY{jj}{kk,catIdx}.Mean{9}(:,xIdx);
            COHR.FlySTD{ww}(jj,kk)      = HeadFree{ww}.FLY{jj}{kk,catIdx}.STD{9}(:,xIdx);
        end
        FREQ.GrandMean{ww}(jj,1) 	= HeadFree{ww}.GRAND{jj,catIdx}.Mean{2}{4};
        COHR.GrandMean{ww}(jj,1)   	= HeadFree{ww}.GRAND{jj,catIdx}.Mean{2}{9}(:,xIdx);
        COHR.GrandSTD{ww}(jj,1)   	= HeadFree{ww}.GRAND{jj,catIdx}.STD{2}{9}(:,xIdx);
    end
end

FIG = figure (figNum); % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 680 500];
FIG.Name = filename;
movegui(FIG,'center')
for ww = 1:nAmp
   FIG.Name = [FIG.Name '_' num2str(Amp(ww))];  
end
hold on

% Fly Stats
for ww = 1:nAmp % amplitudes
    for kk = 1:HeadFree{ww}.N{1,1} % flys
        hold on ; xlim([0 12.5]) ; ylim(1*[0 1])
        h.Fly = errorbar(FREQ.FlyMean{ww}(:,kk),COHR.FlyMean{ww}(:,kk),COHR.FlySTD{ww}(:,kk),'-','LineWidth',1);
        h.Fly.Color(4) = 0.5;
    end
end

% Grand Stats
for ww = 1:nAmp % amplitudes
	hold on ; xlim([0 12.5]) ; ylim(1*[0 1]) ; title([num2str(Amp) , char(176)])
    h.Fly = errorbar(FREQ.GrandMean{ww},COHR.GrandMean{ww},2*COHR.GrandSTD{ww},'-ok','LineWidth',3);
    h.Fly.Color(4) = 0.5;
    ylabel('Coherence')
    xlabel('Frequency (Hz)')
end

figdir = 'H:\DATA\Rigid_Data\FIGURE\';
saveas(FIG,[figdir FIG.Name '.fig']); % save .fig file
% print (FIG,[figdir FIG.Name],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
end