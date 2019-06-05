function [FIG] = MakeFig_Sine_HeadFree_TimeDelay()
%% MakeFig_Sine_HeadFree_TimeDelay:
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
    HeadFree{ww} = load(fullfile(root,FILES{ww}),'TRIAL','U','N');
end

figNum = 1;
filename = 'Sine_HeadFree_TimeDelay';
catIdx = [5 8 7];
xIdx = 1;

nFreq = HeadFree{1}.N{1,3};
Freq = HeadFree{1}.U{1,3}{1};
TD = cell(nFreq,nAmp);
PD = cell(nFreq,nAmp);
ALL = cell(1,nAmp);
for ww = 1:nAmp
    for jj = 1:HeadFree{ww}.N{1,3}
        pp = 1;
        for kk = 1:HeadFree{ww}.N{1,1}
            for ii = 1:size(HeadFree{ww}.TRIAL{kk,jj},1)
                for qq = 1:length(catIdx)
                    TD{jj,ww}(pp,qq) = HeadFree{ww}.TRIAL{kk,jj}{ii,catIdx(qq)}.TimeDiff(:,xIdx);
                    
                    IOFreq = HeadFree{ww}.TRIAL{kk,jj}{ii,catIdx(qq)}.IOFreq;
                    IOPhase = HeadFree{ww}.TRIAL{kk,jj}{ii,catIdx(qq)}.IOBodePhaseDiff(:,xIdx);
                  	PD{jj,ww}(pp,qq) = (IOPhase/(2*pi))*(1/IOFreq);
                end
                pp = pp + 1;
            end
        end
        nTrial = size(TD{jj,ww},1);
        ALL{ww} = [ALL{ww} ; [TD{jj,ww} , PD{jj,ww}, jj*ones(nTrial,1)]];
    end
end

FIG = figure (figNum) ; clf
FIG.Color = 'w';
FIG.Position = [100 100 1200 700];
FIG.Name = filename;
movegui(FIG,'center')

for ww = 1:nAmp
   FIG.Name = [FIG.Name '_' num2str(Amp(ww))];  
end

CC = jet(nFreq);
Labels = {'Ref to Head','Ref to Wings','Head to Wings'};
hold on
for qq = 1:2*length(catIdx)
    ax = subplot(2,length(catIdx),qq); hold on
        ax.FontSize = 12;
        ax.XLabel.FontSize = 14; 
        ax.YLabel.FontSize = ax.XLabel.FontSize;
        
        if qq<=length(catIdx)
            ax.Title.String = Labels{qq};
            ax.XTickLabels = '';
        else
            ax.XLabel.String = 'Frequency (Hz)';
        end

        if qq==1
            ax.YLabel.String = 'CC Time Delay (s)';
        elseif qq==length(catIdx)+1
            ax.YLabel.String = 'Phase Time Delay (s)';
        end
        
        bx = boxplot(ALL{1}(:,qq),ALL{1}(:,end),'Labels',{num2str(Freq)},...
            'Width',0.5,'Symbol','+','Whisker',2);

        h = get(bx(5,:),{'XData','YData'});
        for kk = 1:size(h,1)
           patch(h{kk,1},h{kk,2},CC(kk,:))
        end

        set(findobj(gcf,'tag','Median'), 'Color', 'k','LineWidth',2);
        set(findobj(gcf,'tag','Box'), 'Color', 'k','LineWidth',2);
        set(findobj(gcf,'tag','Upper Whisker'), 'Color', 'k');
        ax.Children = ax.Children([end 1:end-1]);
        ax.YLim = [-0.2 0.2];
end

end