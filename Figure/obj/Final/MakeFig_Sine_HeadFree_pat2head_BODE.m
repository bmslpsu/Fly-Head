function [FIG] = MakeFig_Sine_HeadFree_pat2head_BODE()
%% MakeFig_Sine_HeadFree_pat2head_BODE:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     : figure handle
%

root = 'H:\DATA\Rigid_Data\';

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
    HeadFree{ww} = load(fullfile(root,FILES{ww}),'TRIAL','FLY','GRAND','U','N','T');
end
%%
clearvars -except HeadFree Amp nAmp
figNum = 1;
filename = 'Sine_HeadFree_pat2head_BODE'; % name of figure to save
catIdx = 5; % pat2head
xIdx = 1;
clc
% Fly Stats
FREQ.FlyMean    = cell(nAmp,1);
GAIN.FlyMean  	= cell(nAmp,1);
PHASE.FlyMean  	= cell(nAmp,1);
FREQ.GrandMean	= cell(nAmp,1);
GAIN.GrandMean 	= cell(nAmp,1);
PHASE.GrandMean	= cell(nAmp,1);
GAIN.FlySTD     = cell(nAmp,1);
PHASE.FlySTD  	= cell(nAmp,1);
GAIN.GrandSTD	= cell(nAmp,1);
PHASE.GrandSTD	= cell(nAmp,1);
for ww = 1:nAmp % amplitudes
    for jj = 1:HeadFree{ww}.N{1,3} % frequencies
        for kk = 1:HeadFree{ww}.N{1,1} % flies
            FREQ.FlyMean{ww}(jj,kk)     = HeadFree{ww}.FLY{jj}{kk,catIdx}.Mean{4};
            GAIN.FlyMean{ww}(jj,kk)     = HeadFree{ww}.FLY{jj}{kk,catIdx}.Mean{5}(:,xIdx);
            PHASE.FlyMean{ww}(jj,kk)    = HeadFree{ww}.FLY{jj}{kk,catIdx}.CircMean{6}(:,xIdx);
            GAIN.FlySTD{ww}(jj,kk)      = HeadFree{ww}.FLY{jj}{kk,catIdx}.STD{5}(:,xIdx);
            PHASE.FlySTD{ww}(jj,kk)     = HeadFree{ww}.FLY{jj}{kk,catIdx}.CircSTD{6}(:,xIdx);
            
            if jj==5 && PHASE.FlyMean{ww}(jj,kk)>deg2rad(50)
                PHASE.FlyMean{ww}(jj,kk) = PHASE.FlyMean{ww}(jj,kk) - 1*pi;
            end
            
%             if jj==5 && PHASE.FlyMean{ww}(jj,kk)<deg2rad(-150)
%                 PHASE.FlyMean{ww}(jj,kk) = PHASE.FlyMean{ww}(jj,kk) + 1*pi;
%                 disp('here')
%             end
            
        end
        FREQ.GrandMean{ww}(jj,1) 	= HeadFree{ww}.GRAND{jj,catIdx}.Mean{2}{4};
        GAIN.GrandMean{ww}(jj,1)   	= HeadFree{ww}.GRAND{jj,catIdx}.Mean{2}{5}(:,xIdx);
        PHASE.GrandMean{ww}(jj,1) 	= HeadFree{ww}.GRAND{jj,catIdx}.CircMean{7}{6}(:,xIdx);
        GAIN.GrandSTD{ww}(jj,1)   	= HeadFree{ww}.GRAND{jj,catIdx}.STD{2}{5}(:,xIdx);
        PHASE.GrandSTD{ww}(jj,1) 	= circ_mean(PHASE.FlySTD{ww}(jj,:),[],2);
    end
end

FIG = figure (figNum); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 4 4];
FIG.Name = filename;
movegui(FIG,'center')
for ww = 1:nAmp
   FIG.Name = [FIG.Name '_' num2str(Amp(ww))];  
end
hold on

% % Fly Stats
% for ww = 1:nAmp % amplitudes
%     for kk = 1:HeadFree{ww}.N{1,1} % flys
%       	subplot(2,1,1) ; hold on ; xlim([0 12.5]) ; ylim(1*[0 1])
%             h.Fly = errorbar(FREQ.FlyMean{ww}(:,kk),GAIN.FlyMean{ww}(:,kk),GAIN.FlySTD{ww}(:,kk),'-','LineWidth',1);
%             h.Fly.Color(4) = 0.5;
%         subplot(2,1,2) ; hold on ; xlim([0 12.5]) ; ylim(pi*[-1 1])
%             h.Fly = errorbar(FREQ.FlyMean{ww}(:,kk),rad2deg(PHASE.FlyMean{ww}(:,kk)),...
%                 PHASE.FlySTD{ww}(:,kk),'-','LineWidth',1);
%             h.Fly.Color(4) = 0.5;
%     end
% end

% Grand Stats
for ww = 3 % amplitudes
    ax1 = subplot(2,1,1) ; hold on ; xlim([0 12.5]) ; ylim(1*[0 1]) ; title([num2str(Amp(ww)) , char(176)])
    ax1.FontSize = 8;
    ax1.YLabel.FontSize = 8;
    ax1.XLabel.FontSize = 8;
        h.Fly = errorbar(FREQ.GrandMean{ww},GAIN.GrandMean{ww},2*GAIN.GrandSTD{ww},'-ok','LineWidth',2);
        h.Fly.Color(4) = 0.5;
        ylabel(['Gain (' char(176) '/' char(176) ')'])
    ax2 = subplot(2,1,2) ; hold on ; xlim([0 12.5]) ; ylim([-150 150])
    ax2.FontSize = 8;
	ax2.YLabel.FontSize = 8;
    ax2.XLabel.FontSize = 8;
        h.Fly = errorbar(FREQ.GrandMean{ww},rad2deg(PHASE.GrandMean{ww}),...
            rad2deg(2*PHASE.GrandSTD{ww}),'-ok','LineWidth',2);
        h.Fly.Color(4) = 0.5;
        ylabel('Phase Difference (deg)')
        xlabel('Frequency (Hz)')
end

end