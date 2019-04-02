function [] = MakeFig_Sine_HeadFree_Wing_Pos_BODE(amp,root,figNum)
%% MakeFig_Chirp_HeadFree_Head_Pos_BODE:
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
filename = 'Sine_HeadFree_Wing_Pos_BODE'; % name of figure to save
HeadFree = cell(length(amp),1);
for jj = 1:length(amp)
    HeadFree{jj} = load([root 'Sine_HeadFree_' num2str(3.75*amp(jj)) '_DATA_.mat'],...
        'PAT','WING','HEAD','BODE','CROSS','D','I','N','U');
end

legList = string(num2cell(3.75*amp));

FIG = figure (figNum); % figure handle
FIG.Color = 'w';
set(gcf,'Position',[200 200 (1500/5)*length(amp) 700])
for jj = 1:length(amp)
    for kk = 1:HeadFree{amp(jj)}.N{1,1}
        for ii = 1:HeadFree{amp(jj)}.N{1,3}
            subplot(2,length(amp),jj) ; hold on
            h.med = plot( HeadFree{amp(jj)}.U{1,3}{1}(ii) , HeadFree{amp(jj)}.WING.GAIN{kk}{ii} ,'*','Color',[0.5 0.5 0.5],...
                'LineWidth',1);
            subplot(2,length(amp),jj+length(amp)) ; hold on
            h.med = plot( HeadFree{amp(jj)}.U{1,3}{1}(ii) , HeadFree{amp(jj)}.WING.PhaseDiff{kk}{ii} ,'*','Color',[0.5 0.5 0.5],...
                'LineWidth',1);            
        end
        % MAGNITUDE
        subplot(2,length(amp),jj) ; hold on
            SE = 2*HeadFree{amp(jj)}.WING.FlySTD.GAIN{kk}/1;
            SE = SE.*0;
            h.med = errorbar( HeadFree{amp(jj)}.U{1,3}{1} , HeadFree{amp(jj)}.WING.FlyMed.GAIN{kk} , SE ,'-o',...
                'LineWidth',1);
        subplot(2,length(amp),jj+length(amp)) ; hold on
        % PHASE    
            SE = (2*HeadFree{amp(jj)}.WING.FlySTD.PhaseDiff{kk}/1);
            h.med = errorbar( HeadFree{amp(jj)}.U{1,3}{1} , (HeadFree{amp(jj)}.WING.FlyMed.PhaseDiff{kk}) , SE ,'LineWidth',1);
    end
end

for jj = 1:length(amp)
    % MAGNITUDE
	subplot(2,length(amp),jj) ; hold on
    title([num2str(legList{jj}) '$^{\circ}$'],'Interpreter','latex','FontSize',15)
        SE = 2*HeadFree{amp(jj)}.WING.GrandSTD.GAIN/HeadFree{amp(jj)}.N{1,1};
        h.med = errorbar( HeadFree{amp(jj)}.U{1,3}{1} , HeadFree{amp(jj)}.WING.GrandMed.GAIN , SE ,'k','LineWidth',3);
        
     	xlim([0 12])
        ylim([0 0.3])
        if jj==1
            ylabel('Head Gain(${\circ}/{\circ}$)','Interpreter','latex','FontSize',15)
        end
        if jj~=1
            yticks(0)
            yticklabels('')
        end
        xticks(0)
        xticklabels('')
    % PHASE    
    subplot(2,length(amp),jj+length(amp)) ; hold on
        MED = HeadFree{amp(jj)}.WING.GrandMed.PhaseDiff;
        SE = 3*HeadFree{amp(jj)}.WING.GrandSTD.PhaseDiff/HeadFree{amp(jj)}.N{1,1};
        h.phase = errorbar( HeadFree{amp(jj)}.U{1,3}{1} , MED , SE ,'k','LineWidth',3);
        plot([0 15],[0 0],'Color',[0.5 0.5 0.5],'LineStyle','--')
     	xlim([0 12.5])
        ylim([-4 4])
        if jj==1
            ylabel('Head Phase(rad)','Interpreter','latex','FontSize',15)
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