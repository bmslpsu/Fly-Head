function [] = MakeFig_Sine_HeadFree_Head_Pos_TimeDelay(amp,root,figNum)
%% MakeFig_Chirp_HeadFree_Head_Pos_BODE:
%   INPUTS:
%       root: root directory containing data structure
%       figNum: figure #
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE INPUT %
amp = 2:5;
root = 'H:\EXPERIMENTS\Experiment_Sinusoid\DATA\';
figNum = 1;
%---------------------------------------------------------------------------------------------------------------------------------
filename = 'Sine_HeadFree_Head_Pos_TimeDelay'; % name of figure to save
HeadFree{1} = load([root 'Sine_HeadFree_3.75_DATA_.mat'],'PAT','WING','HEAD','BODE','n','unq');
HeadFree{2} = load([root 'Sine_HeadFree_7.5_DATA_.mat'],'PAT','WING','HEAD','BODE','n','unq');
HeadFree{3} = load([root 'Sine_HeadFree_11.25_DATA_.mat'],'PAT','WING','HEAD','BODE','n','unq');
HeadFree{4} = load([root 'Sine_HeadFree_15_DATA_.mat'],'PAT','WING','HEAD','BODE','n','unq');
HeadFree{5} = load([root 'Sine_HeadFree_18.75_DATA_.mat'],'PAT','WING','HEAD','BODE','n','unq');

legList = (string(num2cell(3.75*amp)));

FIG = figure (figNum); % figure handle
FIG.Color = 'w';
set(gcf,'Position',[200 200 1000 (1000/5)*length(amp)])

for jj = 1:length(amp)
    % Pat to Head TIME DELAY
	subplot(length(amp),3,jj + 2*jj - 2) ; hold on
    ylabel({[num2str(legList{jj}) '$^{\circ}$'],'(s)'},'Interpreter','latex','FontSize',15)
    T = 1./HeadFree{amp(jj)}.unq.Freq; % period
    MED = (HeadFree{amp(jj)}.HEAD.GrandMed.PhaseDiff/(2*pi))' .* T;
    SE = (3*HeadFree{amp(jj)}.HEAD.GrandSTD.PhaseDiff/1/(2*pi))' .* T;
    h.phase = errorbar( HeadFree{amp(jj)}.unq.Freq , MED , SE ,'b','LineWidth',3);
    plot([0 15],[0 0],'Color',[0.5 0.5 0.5],'LineStyle','--')
    ylim([-0.1 0.1])
    xlim([0 12.5])
    xticks([0.1 2:2:12])
    if jj==length(amp)
        xlabel('Frequency (Hz)','Interpreter','latex','FontSize',15)
    elseif jj==1
        title('Pat to Head','Interpreter','latex','FontSize',15)
    end
    
    % Pat to Wing TIME DELAY
	subplot(length(amp),3,jj + 2*jj - 1) ; hold on
    T = 1./HeadFree{amp(jj)}.unq.Freq; % period
    MED = (HeadFree{amp(jj)}.WING.GrandMed.PhaseDiff/(2*pi))' .* T;
    SE = (3*HeadFree{amp(jj)}.WING.GrandSTD.PhaseDiff/1/(2*pi))' .* T;
    h.phase = errorbar( HeadFree{amp(jj)}.unq.Freq , MED , SE ,'r','LineWidth',3);
    plot([0 15],[0 0],'Color',[0.5 0.5 0.5],'LineStyle','--')
    ylim([-0.3 0.3])
    xlim([0 12.5])
    xticks([0.1 2:2:12])
    if jj==length(amp)
        xlabel('Frequency (Hz)','Interpreter','latex','FontSize',15)
    elseif jj==1
        title('Pat to Wing','Interpreter','latex','FontSize',15)
    end
    
    % Head to Wing TIME DELAY
	subplot(length(amp),3,jj + 2*jj - 0) ; hold on
    T = 1./HeadFree{amp(jj)}.unq.Freq; % period
    MED = (HeadFree{amp(jj)}.BODE.GrandMed.head2wing.PhaseDiff/(2*pi))' .* T;
    SE = (3*HeadFree{amp(jj)}.BODE.GrandSTD.head2wing.PhaseDiff/1/(2*pi))' .* T;
    h.phase = errorbar( HeadFree{amp(jj)}.unq.Freq , MED , SE ,'COlor',[0.5 0.1 0.7],'LineWidth',3);
    plot([0 15],[0 0],'Color',[0.5 0.5 0.5],'LineStyle','--')
    ylim([-0.1 0.1])
    xlim([0 12.5])
    xticks([0.1 2:2:12])
    if jj==length(amp)
        xlabel('Frequency (Hz)','Interpreter','latex','FontSize',15)
    elseif jj==1
        title('Head to Wing','Interpreter','latex','FontSize',15)
    end
end

saveas(FIG,[root 'FIGURE\' filename '.fig']); % save .fig file
print(gcf,[root 'FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
disp('Saved to')
disp(root)
end
