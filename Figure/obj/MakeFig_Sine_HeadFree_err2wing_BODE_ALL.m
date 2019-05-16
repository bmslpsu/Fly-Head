function [FIG] = MakeFig_Sine_HeadFree_err2wing_BODE_ALL()
%% MakeFig_Sine_HeadFree_err2wing_BODE_ALL:
%   INPUTS:
%      CLR      :   line color
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'F:\DATA\Rigid_Data\';

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
    HeadFree{ww} = load(fullfile(root,FILES{ww}),'GRAND','U','N');
end
%%
figNum = 1;
filename = 'Sine_HeadFree_err2wing_BODE_ALL'; % name of figure to save
catIdx = 6; % pat2head
xIdx = 1;

% Fly Stats
FREQ.GrandMean	= cell(nAmp,1);
GAIN.GrandMean 	= cell(nAmp,1);
PHASE.GrandMean	= cell(nAmp,1);
GAIN.GrandSTD 	= cell(nAmp,1);
PHASE.GrandSTD	= cell(nAmp,1);
for ww = 1:nAmp % amplitudes
    for jj = 1:HeadFree{ww}.N{1,3} % frequencies
        FREQ.GrandMean{ww}(jj,1) 	= HeadFree{ww}.GRAND{jj,catIdx}.Mean{2}{4};
        GAIN.GrandMean{ww}(jj,1)   	= HeadFree{ww}.GRAND{jj,catIdx}.Mean{2}{5}(:,xIdx);
        PHASE.GrandMean{ww}(jj,1) 	= HeadFree{ww}.GRAND{jj,catIdx}.CircMean{9}{6}(:,xIdx);
        GAIN.GrandSTD{ww}(jj,1)   	= HeadFree{ww}.GRAND{jj,catIdx}.STD{2}{5}(:,xIdx);
        PHASE.GrandSTD{ww}(jj,1) 	= HeadFree{ww}.GRAND{jj,catIdx}.CircSTD{10}{6}(:,xIdx);
        if jj==5 || jj==6
            if PHASE.GrandMean{ww}(jj,1)>0.4
                PHASE.GrandMean{ww}(jj,1) 	= PHASE.GrandMean{ww}(jj,1) - pi;
            end
        end

    end
end

FIG = figure (figNum); clf % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 680 680];
FIG.Name = filename;
movegui(FIG,'center')
hold on

CLR         = jet(nAmp);
legList     = cellstr(num2str(Amp))';
hh          = cell(nAmp,1);
% Grand Stats
for ww = 1:nAmp % amplitudes
    subplot(2,1,1) ; hold on ; xlim([0 12.5]) ; ylim(1*[0 0.4]) ; title('All','interpreter','latex','FontSize',18)
        h.Fly = errorbar(FREQ.GrandMean{ww},GAIN.GrandMean{ww},2*GAIN.GrandSTD{ww},'-o','Color',CLR(ww,:),'LineWidth',2);
        ylabel('Wing Gain ($^{\circ}/^{\circ}$)','interpreter','latex','FontSize',15)
    subplot(2,1,2) ; hold on ; xlim([0 12.5]) ; ylim(pi*[-1 1])
        h.Fly = errorbar(FREQ.GrandMean{ww},PHASE.GrandMean{ww},2*PHASE.GrandSTD{ww},'-o','Color',CLR(ww,:),'LineWidth',2);
        hh{ww} = h.Fly;
        plot([0 12],[0 0],'--k')
        ylabel('Wing Phase (rad)','interpreter','latex','FontSize',15)
        xlabel('Frequency (Hz)','interpreter','latex','FontSize',15)
end
legend(cat(1,hh{:}),legList)

figdir = 'H:\DATA\Rigid_Data\FIGURE\';
saveas(FIG,[figdir FIG.Name '.fig']); % save .fig file
% print (FIG,[figdir FIG.Name],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
end