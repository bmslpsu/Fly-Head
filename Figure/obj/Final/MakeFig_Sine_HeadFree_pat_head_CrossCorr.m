function [] = MakeFig_Sine_HeadFree_pat_head_CrossCorr()
%% MakeFig_Sine_HeadFree_pat_head_CrossCorr:
%   INPUTS:
%       -
%   OUTPUTS:
%       -
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
    HeadFree{ww} = load(fullfile(root,FILES{ww}),'TRIAL','GRAND','U','N');
end

%%
clearvars -except HeadFree Amp nAmp

filename = 'Sine_HeadFree_head2wing_CrossCorr';

hold on
catIdx = 8;
xIdx = 1;

FIG = figure (1); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 6 4];
FIG.Name = filename;
movegui(FIG,'center')

for ww = 1:nAmp % amplitudes 
    pp = 1;
    for jj = 1:HeadFree{ww}.N{1,3} % frequencies
        freq = HeadFree{ww}.U{1,3}{1}(jj);
        ax = subplot(ceil(HeadFree{ww}.N{1,3}/1),1,pp) ; hold on
      	ax.FontSize = 8;
%         title([num2str(freq) ' Hz'],'FontSize',8,'FontWeight','bold')
        ax.YLabel.FontSize = 8;
        ax.XLim = 0.1*[-1 1];
        ax.YLim = [-1 5]*10^(4);
        ax.YTickLabels = '';
        for kk = 1:HeadFree{ww}.N{1,1} % flies
            for ii = 1:size(HeadFree{ww}.TRIAL{kk,jj},1) % trials
                CC = HeadFree{ww}.TRIAL{kk,jj}{ii,catIdx}.CrossCorr(:,xIdx);
                TL = HeadFree{ww}.TRIAL{kk,jj}{ii,catIdx}.TimeLags(:,xIdx);
                maxCC = HeadFree{ww}.TRIAL{kk,jj}{ii,catIdx}.MaxCC(:,xIdx);
                TD = HeadFree{ww}.TRIAL{kk,jj}{ii,catIdx}.TimeDiff(:,xIdx);
                
                h.trial = plot(TL,CC,'LineWidth',0.5);
                h.trial.Color = [0 0 0 1];
                
               	plot([TD TD],[ax.YLim(1) maxCC],'r','LineWidth',0.2)
                plot(TD,maxCC,'.','Color','r','MarkerSize',10)
            end
        end
    	CC = HeadFree{ww}.GRAND{jj,catIdx}.Mean{1}{10}(:,xIdx);
        CC_STD = HeadFree{ww}.GRAND{jj,catIdx}.STD{2}{10}(:,xIdx);
        TL = HeadFree{ww}.GRAND{jj,catIdx}.Mean{1}{11}(:,xIdx);
        PlotPatch(CC ,CC_STD, TL, 2, HeadFree{ww}.N{1,1}, [0.4 0 0.8], [0.4 0.4 0.6], 0.5, 2);
        
        TD_med = HeadFree{ww}.GRAND{jj,catIdx}.Median{2}{13}(:,xIdx);
        loc = [mean(ax.XLim) mean(ax.YLim)];
        h.text = text(loc(1),loc(2),[num2str(1000*TD_med) ' ms'],'Color','b','FontSize',8);
        pp = pp + 1;
    end
end

end