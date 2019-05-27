function [] = MakeFig_Sine_HeadFree_pat_head_CrossCorr()
%% MakeFig_Sine_HeadFree_pat_head_CrossCorr:
%   INPUTS:
%       -
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
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
    HeadFree{ww} = load(fullfile(root,FILES{ww}),'TRIAL','FLY','GRAND','U','N');
end
%%
filename = 'Sine_HeadFree_pat_head_CrossCorr';

hold on
catIdx = 5;
xIdx = 1;
figNum = 1;

FIG = figure (figNum); clf
FIG.Color = 'w';
FIG.Position = [100 100 1100 800];
movegui(FIG,'center')
FIG.Name = filename;
for ww = 1:nAmp
   FIG.Name = [FIG.Name '_' num2str(Amp(ww))];  
end

for ww = 1:nAmp % amplitudes 
    pp = 1;
    for jj = 1:HeadFree{ww}.N{1,3} % frequencies
        freq = HeadFree{ww}.U{1,3}{1}(jj);
     	subplot(HeadFree{ww}.N{1,3},1,pp)  ; hold on
     	title([num2str(freq) ' Hz'],'FontSize',14,'FontWeight','bold')
        
        ax.L = gca;
        ax.L.YColor = [0 0 0];
        ax.L.YLabel.String = 'Cross Corr';
        ax.L.YLabel.FontSize = 14;
        ax.L.FontSize = 12;
%         ax.YLim = 1000*[-1 1];
%         ax.XLim = 20*[0 1];
        if pp==HeadFree{ww}.N{1,3}
            ax.L.XLabel.String = 'Time (s)';
            ax.L.XLabel.FontSize = 14;
        else
            ax.L.XTickLabels = [];
        end
        
        h.patch = PlotPatch(HeadFree{ww}.GRAND{jj,catIdx}.Mean{2}{10}(:,xIdx),HeadFree{ww}.GRAND{jj,catIdx}.STD{2}{10}(:,xIdx),...
            HeadFree{ww}.GRAND{jj,catIdx}.Mean{2}{11}(:,xIdx),2,HeadFree{ww}.N{1,1},'b',[0.4 0.4 0.6],0.5,3);
        
       pp = pp + 1;
    end
end

TimeDiff = cell(nAmp,HeadFree{ww}.N{1,3});
for ww = 1:nAmp % amplitudes 
    for jj = 1:HeadFree{ww}.N{1,3} % frequencies
        for kk = 1:HeadFree{ww}.N{1,1} % flies
            for ii = 1:size(HeadFree{ww}.TRIAL{kk,jj},1) % trials
                TimeDiff{ww,jj}(end+1,1) = HeadFree{ww}.TRIAL{kk,jj}{ii,catIdx}.TimeDiff(:,xIdx);
            end
        end
    end
end

BOX = figure (figNum+1); clf
BOX.Color = 'w';
BOX.Position = [100 100 1100 800];
movegui(BOX,'center')


end