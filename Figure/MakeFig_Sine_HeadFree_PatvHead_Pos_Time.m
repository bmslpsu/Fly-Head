function [] = MakeFig_Sine_HeadFree_PatvHead_Pos_Time(amp,root,figNum)
%% MakeFig_Sine_HeadFree_HeadvWing_Pos_Time: time domain head vs wing for all amplitudes & frequencies
%   INPUTS:
%       amp     : amplitudes to plot
%       root    : root directory containing data structure
%       figNum  : figure #
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE INPUT %
root = 'H:\EXPERIMENTS\Experiment_Sinusoid\DATA\';
figNum = 1;
amp = 1:5;
%---------------------------------------------------------------------------------------------------------------------------------
filename = 'Sine_HeadFree_PatvHead_Pos_Time'; % name of figure to save
HeadFree{1} = load([root 'Sine_HeadFree_3.75_DATA_.mat'],'PAT','WING','HEAD','BODE','n','unq');
HeadFree{2} = load([root 'Sine_HeadFree_7.5_DATA_.mat'],'PAT','WING','HEAD','BODE','n','unq');
HeadFree{3} = load([root 'Sine_HeadFree_11.25_DATA_.mat'],'PAT','WING','HEAD','BODE','n','unq');
HeadFree{4} = load([root 'Sine_HeadFree_15_DATA_.mat'],'PAT','WING','HEAD','BODE','n','unq');
HeadFree{5} = load([root 'Sine_HeadFree_18.75_DATA_.mat'],'PAT','WING','HEAD','BODE','n','unq');

close all
FIG = figure (figNum); % figure handle
FIG.Color = 'w';
FIG.Position = [100 0 1400 900];
Fs = 200;
T = 10;
R = nan(length(amp),HeadFree{amp(1)}.n.Freq);
for kk = 1:length(amp)
    for jj = 1:HeadFree{amp(kk)}.N{1,1}
        subIdx = jj + HeadFree{amp(kk)}.N{1,1}*(kk-1);
        subplot(length(amp),HeadFree{amp(kk)}.N{1,1},subIdx) ; hold on
        if subIdx<=HeadFree{amp(kk)}.N{1,1}
        	title([num2str(HeadFree{amp(kk)}.U{1,3}{1}(jj)) ' Hz'],'Interpreter','latex','FontSize',12)
        end
        if subIdx>=(HeadFree{amp(kk)}.N{1,1})*(length(amp)-1)
            xlabel('Stimulus ($^{\circ}$)','Interpreter','latex','FontSize',15)
            xticks(-20:10:20)
        else
            xticks(0)
            xticklabels('')
        end
        if (subIdx-1)/HeadFree{amp(kk)}.N{1,1} == round((subIdx-1)/HeadFree{amp(kk)}.N{1,1})
        	ylabel({[num2str(3.75*amp(kk)) '$^{\circ}$'],'Head ($^{\circ}$)'},'Interpreter','latex','FontSize',15)
            yticks(-20:10:20)
        else
            yticks(0)
            yticklabels('')
        end
        xData = HeadFree{amp(kk)}.PAT.ALL.Pos{jj}(1*Fs:T*Fs,:);
        yData = HeadFree{amp(kk)}.HEAD.ALL.Pos{jj}(1*Fs:T*Fs,:);
        [h] = scatplot(xData(:),yData(:));
        axis([-20 20 -10 10])
        delete(h.cb)

        % Calculate linear best fit
        [r,m,b] = regression(xData(:),yData(:),'one');
        R(kk,jj) = r;
        text(-19,9,['r =' num2str(r)]);
        xFit = linspace(-20,20,4000);
        yFit = m*xFit + b;
        plot(xFit,yFit,'r','LineWidth',5)

    end
end
% saveas(FIG,[root 'FIGURE\' filename '.fig']); % save .fig file
% print(gcf,[root 'FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
% disp('Saved to')
% disp(root)
end