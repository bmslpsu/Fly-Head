function [] = MakeFig_Sine_HeadFree_Pat2Head_Pos_TimeDelay(amp,root,figNum)
%% MakeFig_Sine_HeadFree_Head_Pos_Time: time domain sinusoid plot for all amplitudes & frequencies
%   INPUTS:
%       amp     : amplitudes to plot
%       root    : root directory containing data structure
%       figNum  : figure #
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE INPUT %
root = 'F:\EXPERIMENTS\Experiment_Sinusoid\DATA\';
figNum = 1;
amp = 4;
%---------------------------------------------------------------------------------------------------------------------------------
filename = 'Sine_HeadFree_Pat2Head_Pos_TimeDelay'; % name of figure to save
HeadFree{1} = load([root 'Sine_HeadFree_3.75_DATA_.mat'],'PAT','WING','HEAD','BODE','CROSS','n','unq');
HeadFree{2} = load([root 'Sine_HeadFree_7.5_DATA_.mat'],'PAT','WING','HEAD','BODE','CROSS','n','unq');
HeadFree{3} = load([root 'Sine_HeadFree_11.25_DATA_.mat'],'PAT','WING','HEAD','BODE','CROSS','n','unq');
HeadFree{4} = load([root 'Sine_HeadFree_15_DATA_.mat'],'PAT','WING','HEAD','BODE','CROSS','n','unq');
HeadFree{5} = load([root 'Sine_HeadFree_18.75_DATA_.mat'],'PAT','WING','HEAD','BODE','CROSS','n','unq');
%%

% FIG = figure (figNum); % figure handle
% FIG.Color = 'w';
% FIG.Position = [100 100 1200 850];

G = [];
catg = 1:HeadFree{1}.n.Freq;
for jj = catg
    [rr,~] = find( HeadFree{4}.CROSS.ALL(:,3) == HeadFree{1}.unq.Freq(jj) );
    G(rr,1) = jj;
end

SS = [6,7,8]; % box plot data column indicies
YY = { 'Pat to Head (ms)', 'Pat to Wing (ms)' , 'Head to Wing (ms)' }; % ylabels
    
CC = {[0 0 0.5],[0 0 0],[0 0.7 0 ],[0.7 0 0],[0.5 0.5 0.5],[0.1 0.5 0.5],[0.5 0.1 0.5],[0.9 0.7 0.7],[0.3 0.3 0.9],...
    [0.6 0.1 0.8],[0 0.8 0.8]};
FF = cell(length(SS),1);
AX = cell(length(SS),1);
for ww = 1:length(SS)
    FF{ww} = figure (400+ww) ; clf ; hold on
    FF{ww}.Position = [-100+(100*ww) , 500-(20*ww) , 700 , 500];
    FF{ww}.Color = 'w';
    AX{ww} = gca;
    axis tight
        bx = boxplot(HeadFree{4}.CROSS.ALL(:,SS(ww)),G,'Labels',{num2str(HeadFree{4}.unq.Freq(catg,:))},...
            'Width',0.5,'Symbol','+','Whisker',2);
        ax = gca;
        xlabel('Frequency (Hz)')
        ylabel(YY{ww})
        set(gca,'FontSize',20)
        h = get(bx(5,:),{'XData','YData'});
        for k=1:size(h,1)
           patch(h{k,1},h{k,2},CC{ww});
        end
        set(findobj(gcf,'tag','Median'), 'Color', 'w');
        set(findobj(gcf,'tag','Box'), 'Color', 'k');
        set(findobj(gcf,'tag','Upper Whisker'), 'Color', 'k');
        ax.Children = ax.Children([end 1:end-1]);
end












% saveas(FIG,[root 'FIGURE\' filename '.fig']); % save .fig file
% print(gcf,[root 'FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
% disp('Saved to')
% disp(root)
end