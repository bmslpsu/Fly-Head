function [] = MakeFig_Sine_HeadFree_Head_Pos_BODE(amp,root,figNum)
%% MakeFig_Chirp_HeadFree_Head_Pos_BODE:
%   INPUTS:
%       root: root directory containing data structure
%       figNum: figure #
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE INPUT %
% amp = 2:5;
% root = 'H:\EXPERIMENTS\Experiment_Sinusoid\DATA\';
% figNum = 1;
%---------------------------------------------------------------------------------------------------------------------------------
filename = 'Sine_HeadFree_Head_Pos_BODE'; % name of figure to save
HeadFree{1} = load([root 'Sine_HeadFree_3.75_DATA_.mat'],'PAT','WING','HEAD','BODE','n','unq');
HeadFree{2} = load([root 'Sine_HeadFree_7.5_DATA_.mat'],'PAT','WING','HEAD','BODE','n','unq');
HeadFree{3} = load([root 'Sine_HeadFree_11.25_DATA_.mat'],'PAT','WING','HEAD','BODE','n','unq');
HeadFree{4} = load([root 'Sine_HeadFree_15_DATA_.mat'],'PAT','WING','HEAD','BODE','n','unq');
HeadFree{5} = load([root 'Sine_HeadFree_18.75_DATA_.mat'],'PAT','WING','HEAD','BODE','n','unq');

legList = (string(num2cell(3.75*amp)));

FIG = figure (figNum); % figure handle
FIG.Color = 'w';
set(gcf,'Position',[200 200 (1500/5)*length(amp) 700])
for jj = 1:length(amp)
    for kk = 1:HeadFree{amp(jj)}.n.Fly
        for ii = 1:HeadFree{amp(jj)}.n.Freq
            subplot(2,length(amp),jj) ; hold on
            h.med = plot( HeadFree{amp(jj)}.unq.Freq(ii) , HeadFree{amp(jj)}.HEAD.GAIN{kk}{ii} ,'*','Color',[0.5 0.5 0.5],...
                'LineWidth',1);
            subplot(2,length(amp),jj+length(amp)) ; hold on
            h.med = plot( HeadFree{amp(jj)}.unq.Freq(ii) , HeadFree{amp(jj)}.HEAD.PhaseDiff{kk}{ii} ,'*','Color',[0.5 0.5 0.5],...
                'LineWidth',1);            
        end
        % MAGNITUDE
        subplot(2,length(amp),jj) ; hold on
            SE = 2*HeadFree{amp(jj)}.HEAD.FlySTD.GAIN{kk}/1;
            SE = SE.*0;
            h.med = errorbar( HeadFree{amp(jj)}.unq.Freq , HeadFree{amp(jj)}.HEAD.FlyMed.GAIN{kk} , SE ,'-o',...
                'LineWidth',1);
        subplot(2,length(amp),jj+length(amp)) ; hold on
        % PHASE    
            SE = (2*HeadFree{amp(jj)}.HEAD.FlySTD.PhaseDiff{kk}/1);
            h.med = errorbar( HeadFree{amp(jj)}.unq.Freq , (HeadFree{amp(jj)}.HEAD.FlyMed.PhaseDiff{kk}) , SE ,'LineWidth',1);
    end
end

for jj = 1:length(amp)
    % MAGNITUDE
	subplot(2,length(amp),jj) ; hold on
    title([num2str(legList{jj}) '$^{\circ}$'],'Interpreter','latex','FontSize',15)
        SE = 2*HeadFree{amp(jj)}.HEAD.GrandSTD.GAIN/HeadFree{amp(jj)}.n.Fly;
        h.med = errorbar( HeadFree{amp(jj)}.unq.Freq , HeadFree{amp(jj)}.HEAD.GrandMed.GAIN , SE ,'k','LineWidth',3);
        
     	xlim([0 12])
        ylim([0 1])
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
        MED = HeadFree{amp(jj)}.HEAD.GrandMed.PhaseDiff;
        SE = 3*HeadFree{amp(jj)}.HEAD.GrandSTD.PhaseDiff/HeadFree{amp(jj)}.n.Fly;;
        h.phase = errorbar( HeadFree{amp(jj)}.unq.Freq , MED , SE ,'k','LineWidth',3);
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

% saveas(FIG,[root 'FIGURE\' filename '.fig']); % save .fig file
% print(gcf,[root 'FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
% disp('Saved to')
% disp(root)
end
