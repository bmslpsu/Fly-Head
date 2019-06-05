function [] = MakeFig_Sine_HeadFree_phase_trial()
%% MakeFig_Sine_HeadFree_phase_trial:
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
%%
filename = 'Sine_HeadFree_phase_trial';

hold on
catIdx = 5;
xIdx = 1;
figNum = 1;

FIG = figure (figNum); clf
FIG.Color = 'w';
FIG.Position = [100 100 800 500];
movegui(FIG,'center')
FIG.Name = filename;
for ww = 1:nAmp
   FIG.Name = [FIG.Name '_' num2str(Amp(ww))];  
end

CC = jet(HeadFree{ww}.N{1,1});

freqPlot = 2;
nPlot = length(freqPlot);
for ww = 1:nAmp
    pp = 1;
    for jj = freqPlot
        for kk = 1:HeadFree{ww}.N{1,1}
            nTrial = size(HeadFree{ww}.TRIAL{kk,jj},1);
            trial = 1:nTrial;
            Phase = nan(nTrial,1);
            for ii = 1:nTrial
                freq = HeadFree{ww}.U{1,3}{1}(jj);
                ax = subplot(nPlot,1,pp) ; hold on
                ax.FontSize         = 12;
                ax.Title.String     = [num2str(freq) ' Hz'];
                ax.Title.FontSize   = 14;
                ax.Title.FontWeight = 'bold';
                ax.YLim             = 120*[0 1];
                ax.YLabel.String    = ['Phase (' char(176) ')'];
                ax.YLabel.FontSize  = 14;
                ax.XLim             = [0.8 5.2];
               	ax.XLabel.String    = 'Trial #';
                ax.XLabel.FontSize  = ax.YLabel.FontSize;
                ax.XTick            = 1:5;
                
                if pp~=nPlot
                    ax.XTickLabels = [];
                end
                
            	Phase(ii) = rad2deg(HeadFree{ww}.TRIAL{kk,jj}{ii,catIdx}.IOBodePhaseDiff(:,xIdx));
%                 Phase(Phase>120) = Phase(Phase>120) - 360;
            end
            plot(trial,Phase,'o-','Color',CC(kk,:),'LineWidth',2)
        end
        pp = pp + 1;
    end
%     plot([1 5],[0 0],'--k')
end
end