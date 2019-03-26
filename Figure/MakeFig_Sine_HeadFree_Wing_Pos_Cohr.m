function [] = MakeFig_Sine_HeadFree_Wing_Pos_Cohr(amp,root,figNum)
%% MakeFig_Sine_HeadFree_Wing_Pos_Cohr: wing coherence
%   INPUTS:
%       amp     : amplitude indicies to plot
%       root    : root directory containing data structure
%       figNum  : figure #
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE INPUT %
% amp = 1:5;
% root = 'H:\EXPERIMENTS\Experiment_Sinusoid\DATA\';
% figNum = 1;
%---------------------------------------------------------------------------------------------------------------------------------
filename = 'Sine_HeadFree_Wing_Pos_Cohr'; % name of figure to save
HeadFree = cell(length(amp),1);
for jj = 1:length(amp)
    HeadFree{jj} = load([root 'Sine_HeadFree_' num2str(3.75*amp(jj)) '_DATA_.mat'],...
        'PAT','WING','HEAD','BODE','CROSS','D','I','N','U');
end

legList = string(num2cell(3.75*amp));

FIG = figure (figNum); % figure handle
FIG.Color = 'w';
set(gcf,'Position',[200 200 (1500/5)*length(amp) 700/2])
for jj = 1:length(amp)
    for kk = 1:HeadFree{amp(jj)}.N{1,1}
        for ii = 1:HeadFree{amp(jj)}.N{1,3}
            subplot(1,length(amp),jj) ; hold on
            h.med = plot( HeadFree{amp(jj)}.U{1,3}{1}(ii) , HeadFree{amp(jj)}.WING.COHR.Mag{kk}{ii} ,'*','Color',[0.5 0.5 0.5],...
                'LineWidth',1);      
        end
        subplot(1,length(amp),jj) ; hold on
        SE = 2*HeadFree{amp(jj)}.WING.FlySTD.COHR.Mag{kk}/1;
        SE = SE.*0;
        h.med = errorbar( HeadFree{amp(jj)}.U{1,3}{1} , HeadFree{amp(jj)}.WING.FlyMed.COHR.Mag{kk} , SE ,'-o',...
            'LineWidth',1);
    end
end

for jj = 1:length(amp)
    subplot(1,length(amp),jj) ; hold on
    title([num2str(legList{jj}) '$^{\circ}$'],'Interpreter','latex','FontSize',15)
        SE = 2*HeadFree{amp(jj)}.WING.GrandSTD.COHR.Mag;
        h.med = errorbar( HeadFree{amp(jj)}.U{1,3}{1} , HeadFree{amp(jj)}.WING.GrandMed.COHR.Mag , SE ,'k','LineWidth',3);
        
     	xlim([0 12])
        ylim([0 1])
        if jj==1
            ylabel('Head Coherence','Interpreter','latex','FontSize',15)
        end
        if jj~=1
            yticks(0)
            yticklabels('')
        end
        xticks([0.1 2:2:12])
        xlabel('Frequency (Hz)','Interpreter','latex','FontSize',15)
end

saveas(FIG,[root 'FIGURE\' filename '.fig']); % save .fig file
print(gcf,[root 'FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
disp('Saved to')
disp(root)
end
