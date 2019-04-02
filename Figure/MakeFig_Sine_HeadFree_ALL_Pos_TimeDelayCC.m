function [] = MakeFig_Sine_HeadFree_ALL_Pos_TimeDelayCC(amp,root,figNum)
%% MakeFig_Sine_HeadFree_Head_Pos_Time: time domain sinusoid plot for all amplitudes & frequencies
%   INPUTS:
%       amp     : amplitudes to plot
%       root    : root directory containing data structure
%       figNum  : figure #
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE INPUT %
% root = 'H:\EXPERIMENTS\Experiment_Sinusoid\DATA\';
% figNum = 1;
% amp = 2:5;
%---------------------------------------------------------------------------------------------------------------------------------
filename = 'Sine_HeadFree_ALL_Pos_TimeDelayCC'; % name of figure to save
HeadFree{1} = load([root 'Sine_HeadFree_3.75_DATA_.mat'],'PAT','WING','HEAD','BODE','CROSS','n','unq');
HeadFree{2} = load([root 'Sine_HeadFree_7.5_DATA_.mat'],'PAT','WING','HEAD','BODE','CROSS','n','unq');
HeadFree{3} = load([root 'Sine_HeadFree_11.25_DATA_.mat'],'PAT','WING','HEAD','BODE','CROSS','n','unq');
HeadFree{4} = load([root 'Sine_HeadFree_15_DATA_.mat'],'PAT','WING','HEAD','BODE','CROSS','n','unq');
HeadFree{5} = load([root 'Sine_HeadFree_18.75_DATA_.mat'],'PAT','WING','HEAD','BODE','CROSS','n','unq');

FIG = figure (figNum); % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1200 (1000/5)*length(amp)];
SS = [6,7,8]; % box plot data column indicies
YY = { 'Pat to Head', 'Pat to Wing' , 'Head to Wing' }; % ylabels
CC = {[0 0 1],[1 0 0],[0.5 0.1 0.5 ]};
for kk = 1:length(amp)
    G = nan(length(HeadFree{amp(kk)}.CROSS.ALL),1);
    catg = 1:HeadFree{amp(kk)}.n.Freq;
    for jj = catg
        [rr,~] = find( HeadFree{amp(kk)}.CROSS.ALL(:,3) == HeadFree{amp(kk)}.unq.Freq(jj) );
        G(rr,1) = jj;
    end
    
    for ww = 1:length(SS)
        subIdx = ww + length(SS)*(kk-1);
        subplot(length(amp),length(SS),subIdx) ; hold on ; axis tight
        bx = boxplot(HeadFree{amp(kk)}.CROSS.ALL(:,SS(ww)),G,'Labels',{num2str(HeadFree{amp(kk)}.unq.Freq(catg,:))},...
            'Width',0.5,'Symbol','+','Whisker',2);
        ax = gca;
        ax.FontSize = 12;
        
        if subIdx/length(SS)==round(subIdx/length(SS))
            ylim([-0.2 0.2])
        else
            ylim([-2 2])
        end

        if subIdx<=length(SS)
            title(YY{ww})
        end
        
        if subIdx>length(SS)*(length(amp)-1)
            xlabel('Frequency (Hz)')
            vec_pos = get(get(gca, 'XLabel'), 'Position');
            set(get(gca, 'XLabel'), 'Position', vec_pos + [0 -0.01 0]);
        else
            xticks(0)
            xticklabels('')
        end
        
        if (subIdx-1)/length(SS)==round((subIdx-1)/length(SS))
            ylabel({[num2str(3.75*amp(kk)) char(176) ],'s'})
        else
%             yticks(0)
%             yticklabels('')
        end
        h = get(bx(5,:),{'XData','YData'});
        for k=1:size(h,1)
           patch(h{k,1},h{k,2},CC{ww});
        end
        set(findobj(gcf,'tag','Median'), 'Color', 'w');
        set(findobj(gcf,'tag','Box'), 'Color', 'k');
        set(findobj(gcf,'tag','Upper Whisker'), 'Color', 'k');
        ax.Children = ax.Children([end 1:end-1]);
    end
end

saveas(FIG,[root 'FIGURE\' filename '.fig']); % save .fig file
print(gcf,[root 'FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
disp('Saved to')
disp(root)
end