function [FIG] = MakeFig_Sine_HeadFree_pat2head_COHR_ALL()
%% MakeFig_Sine_HeadFree_pat2head_COHR_ALL:
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

figNum = 1;
filename = 'Sine_HeadFree_pat2head_COHR_ALL'; % name of figure to save
catIdx = 5; % pat2head
xIdx = 1;

% Fly Stats
FREQ.GrandMean	= cell(nAmp,1);
COHR.GrandMean 	= cell(nAmp,1);
COHR.GrandSTD 	= cell(nAmp,1);
for ww = 1:nAmp % amplitudes
    for jj = 1:HeadFree{ww}.N{1,3} % frequencies
        FREQ.GrandMean{ww}(jj,1) 	= HeadFree{ww}.GRAND{jj,catIdx}.Mean{2}{4};
        COHR.GrandMean{ww}(jj,1)   	= HeadFree{ww}.GRAND{jj,catIdx}.Mean{2}{9}(:,xIdx);
        COHR.GrandSTD{ww}(jj,1)   	= HeadFree{ww}.GRAND{jj,catIdx}.STD{2}{9}(:,xIdx);

    end
end

FIG = figure (figNum); clf % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 680 500];
FIG.Name = filename;
movegui(FIG,'center')
hold on

CLR         = jet(nAmp);
legList     = cellstr(num2str(Amp))';
hh          = cell(nAmp,1);
% Grand Stats
for ww = 1:nAmp % amplitudes
    hold on ; xlim([0 12.5]) ; ylim(1*[0 1]) ; title('All','interpreter','latex','FontSize',18)
    h.Fly = errorbar(FREQ.GrandMean{ww},COHR.GrandMean{ww},2*COHR.GrandSTD{ww},'-o','Color',CLR(ww,:),'LineWidth',2);
    ylabel('Coherence','interpreter','latex','FontSize',15)
    xlabel('Frequency (Hz)','interpreter','latex','FontSize',15)
end
legend(cat(1,hh{:}),legList)

figdir = 'F:\DATA\Rigid_Data\FIGURE\';
saveas(FIG,[figdir FIG.Name '.fig']); % save .fig file
% print (FIG,[figdir FIG.Name],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
end