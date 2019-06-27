function [FIG] = MakeFig_Sine_HeadFree_CC()
%% MakeFig_Sine_HeadFree_CC:
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
    HeadFree{ww} = load(fullfile(root,FILES{ww}),'GRAND','U','N');
end
%%
figNum = 1;
filename = 'MakeFig_Sine_HeadFree_CC';
catIdx = [5 8 7];
xIdx = 1;

nFreq = HeadFree{1}.N{1,3};
Freq = HeadFree{1}.U{1,3}{1};
TL = cell(nAmp,nFreq);
CC = cell(nAmp,nFreq);
TD = nan(nAmp,nFreq);
MC = nan(nAmp,nFreq);
for ww = 1:nAmp
    for jj = 1:HeadFree{ww}.N{1,3}
        TL{ww,jj} = HeadFree{ww}.GRAND{jj,catIdx(1)}.Mean{1}{11}(:,xIdx);
        CC{ww,jj} = HeadFree{ww}.GRAND{jj,catIdx(1)}.Mean{1}{10}(:,xIdx);
        [~,idx] = max(abs(CC{ww,jj}));
        MC(ww,jj) = CC{ww,jj}(idx);
        TD(ww,jj) = TL{ww,jj}(idx);        
    end
end


FIG = figure (figNum) ; clf
FIG.Color = 'w';
FIG.Position = [100 100 900 500];
FIG.Name = filename;
movegui(FIG,'center')

for ww = 1:nAmp
   FIG.Name = [FIG.Name '_' num2str(Amp(ww))];  
end

cList = jet(nFreq);
Labels = {'Ref to Head','Ref to Wings','Head to Wings'};
hold on

for ww = 1
    pp = 1;
    for jj = 1:HeadFree{ww}.N{1,3}
        ax = subplot(ceil(nFreq/3),3,pp); hold on
        ax.YLabel.String = [num2str(HeadFree{ww}.U{1,3}{1}(jj)) ' Hz'];
        ax.YLabel.FontSize = 14;
        ax.XLim = 0.1*[-1 1];
%         ax.XTick = ax.XLim(1):0.025:ax.XLim(end);
        ax.YTick = [];
        
        if pp>nFreq-3
            ax.XLabel.String = 'Time Difference (s)';
            ax.XLabel.FontSize = ax.YLabel.FontSize;
        end

        plot(TL{ww,jj},CC{ww,jj},'Color','k','LineWidth',2)
        plot(TD(ww,jj),MC(ww,jj),'o','MarkerSize',8)
        plot(TD(ww,jj)*ones(1,2),[-100,MC(ww,jj)],'--r')
        
        text(TD(ww,jj) - 0.05 , MC(ww,jj)-5000 , ['Time Delay = ' num2str(round(1000*TD(ww,jj),1)) 'ms'])
        
        pp = pp + 1;
    end
end

end