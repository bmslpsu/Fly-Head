function [FIG] = MakeFig_Sine_HeadFree_pat2head_BODE_ALL(CLR,saveFig)
%% MakeFig_Sine_HeadFree_pat2head_BODE_ALL:
%   INPUTS:
%      CLR      :   line color
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';
figNum = 1;
filename = 'Sine_HeadFree_pat2head_BODE_ALL'; % name of figure to save
catIdx = 5; % pat2head
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
%%
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
        for kk = 1:HeadFree{ww}.N{1,1} % flys
            FREQ.FlyMean{ww}(jj,kk)     = HeadFree{ww}.FLY{jj}{kk,catIdx}.Mean{4};
            GAIN.FlyMean{ww}(jj,kk)     = HeadFree{ww}.FLY{jj}{kk,catIdx}.Mean{5}(:,xIdx);
            PHASE.FlyMean{ww}(jj,kk)    = HeadFree{ww}.FLY{jj}{kk,catIdx}.CircMean{6}(:,xIdx);
            
            if jj==5
                if PHASE.FlyMean{ww}(5,kk)>0.4
                    PHASE.FlyMean{ww}(5,kk) = PHASE.FlyMean{ww}(5,kk) - pi;
                end
            end
            
            GAIN.FlySTD{ww}(jj,kk)      = HeadFree{ww}.FLY{jj}{kk,catIdx}.STD{5}(:,xIdx);
            PHASE.FlySTD{ww}(jj,kk)     = HeadFree{ww}.FLY{jj}{kk,catIdx}.CircSTD{6}(:,xIdx);
        end
        FREQ.GrandMean{ww}(jj,1) 	= HeadFree{ww}.GRAND{jj,catIdx}.Mean{2}{4};
        GAIN.GrandMean{ww}(jj,1)   	= HeadFree{ww}.GRAND{jj,catIdx}.Mean{2}{5}(:,xIdx);
        PHASE.GrandMean{ww}(jj,1) 	= HeadFree{ww}.GRAND{jj,catIdx}.CircMean{9}{6}(:,xIdx);
        GAIN.GrandSTD{ww}(jj,1)   	= HeadFree{ww}.GRAND{jj,catIdx}.STD{2}{5}(:,xIdx);
        PHASE.GrandSTD{ww}(jj,1) 	= HeadFree{ww}.GRAND{jj,catIdx}.CircSTD{9}{6}(:,xIdx);
        
        if jj==5
            if PHASE.GrandMean{ww}(5,1)>0.4
                PHASE.GrandMean{ww}(5,1) = PHASE.GrandMean{ww}(5,:) - pi;
                PHASE.GrandSTD{ww}(5,1)  = std(PHASE.FlyMean{ww}(5,:));
            end
        end
    end
end

FIG = figure (figNum); % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 680 680];
FIG.Name = filename;
movegui(FIG,'center')
hold on

% Fly Stats
for ww = 1:nAmp % amplitudes
    for kk = 1:HeadFree{ww}.N{1,1} % flys
      	subplot(2,1,1) ; hold on ; xlim([0 12.5]) ; ylim(1*[0 1])
            h.Fly = errorbar(FREQ.FlyMean{ww}(:,kk),GAIN.FlyMean{ww}(:,kk),GAIN.FlySTD{ww}(:,kk),...
                '-','Color',[0.5 0.5 0.5],'LineWidth',1);
            h.Fly.Color(4) = 0.5;
        subplot(2,1,2) ; hold on ; xlim([0 12.5]) ; ylim(pi*[-1 1])
            h.Fly = errorbar(FREQ.FlyMean{ww}(:,kk),PHASE.FlyMean{ww}(:,kk),PHASE.FlySTD{ww}(:,kk),...
                '-','Color',[0.5 0.5 0.5],'LineWidth',1);
            h.Fly.Color(4) = 0.5;
    end
end

% Grand Stats
for ww = 1:nAmp % amplitudes
    subplot(2,1,1) ; hold on ; xlim([0 12.5]) ; ylim(1*[0 1]) ; title('All')
        h.Fly = errorbar(FREQ.GrandMean{ww},GAIN.GrandMean{ww},2*GAIN.GrandSTD{ww},'-o','Color',CLR,'LineWidth',3);
        legend(h.Fly,[num2str(Amp(ww)) char(176)])
        
        ylabel(['Gain (' char(176) '/' char(176) ')'])
    subplot(2,1,2) ; hold on ; xlim([0 12.5]) ; ylim(pi*[-1 1])
        h.Fly = errorbar(FREQ.GrandMean{ww},PHASE.GrandMean{ww},2*PHASE.GrandSTD{ww},'-o','Color',CLR,'LineWidth',3);
        ylabel('Phase Difference (rad)')
        xlabel('Frequency (Hz)')
end

if saveFig
    figdir = 'H:\DATA\Rigid_Data\FIGURE\';
    saveas(FIG,[figdir FIG.Name '.fig']); % save .fig file
    % print (FIG,[figdir FIG.Name],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
end
end