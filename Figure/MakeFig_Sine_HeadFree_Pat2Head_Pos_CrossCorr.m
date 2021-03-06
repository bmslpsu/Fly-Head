function [] = MakeFig_Sine_HeadFree_Pat2Head_Pos_CrossCorr(amp,root,figNum)
%% MakeFig_Sine_HeadFree_Head_Pos_Time: time domain sinusoid plot for all amplitudes & frequencies
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
filename = 'Sine_HeadFree_Pat2Head_Pos_CrossCorr'; % name of figure to save
HeadFree{1} = load([root 'Sine_HeadFree_3.75_DATA_.mat'],'PAT','WING','HEAD','BODE','CROSS','n','unq');
HeadFree{2} = load([root 'Sine_HeadFree_7.5_DATA_.mat'],'PAT','WING','HEAD','BODE','CROSS','n','unq');
HeadFree{3} = load([root 'Sine_HeadFree_11.25_DATA_.mat'],'PAT','WING','HEAD','BODE','CROSS','n','unq');
HeadFree{4} = load([root 'Sine_HeadFree_15_DATA_.mat'],'PAT','WING','HEAD','BODE','CROSS','n','unq');
HeadFree{5} = load([root 'Sine_HeadFree_18.75_DATA_.mat'],'PAT','WING','HEAD','BODE','CROSS','n','unq');
%%
CList = {'k','b','r','g','c'}';
legList = (string(num2cell(3.75*amp)));

FIG = figure (figNum); % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1200 850];
h.med = cell(HeadFree{amp(1)}.n.Freq,length(amp));
for mm = length(amp):-1:1
%     for kk = 1:HeadFree{amp(mm)}.n.Fly
%         for jj = 1:HeadFree{amp(mm)}.n.Freq
%             subplot(HeadFree{amp(mm)}.n.Freq,1,jj) ; hold on
%             h.Trial = plot(HeadFree{amp(mm)}.CROSS.pat2head.lag{kk}{jj},HeadFree{amp(mm)}.CROSS.pat2head.cc{kk}{jj},'Color',[0.5 0.5 0.5],'LineWidth',1); % trials
%             for ii = 1:length(h.Trial)
%                 h.Trial(ii).Color(4) = 0.25;
%             end
%             h.Fly = plot(HeadFree{amp(mm)}.CROSS.FlyMed.pat2head.lag{kk}(:,jj),HeadFree{amp(mm)}.CROSS.FlyMed.pat2head.cc{kk}(:,jj),'LineWidth',2); % flys
%             h.Fly.Color(4) = 0.5;
%         end
%     end
    for jj = 1:HeadFree{amp(mm)}.n.Freq
        subplot(HeadFree{amp(mm)}.n.Freq,1,jj) ; hold on
        title([num2str(HeadFree{amp(mm)}.unq.Freq(jj)) ' Hz'],'Interpreter','latex','FontSize',15)
            [h.patch,h.med{jj,mm}] = PlotPatch(HeadFree{amp(mm)}.CROSS.GrandMed.pat2head.cc(:,jj),HeadFree{amp(mm)}.CROSS.GrandSTD.pat2head.cc(:,jj),...
                HeadFree{amp(mm)}.CROSS.GrandMed.pat2head.lag(:,jj),2,HeadFree{amp(mm)}.n.Fly,CList{amp(mm)},[0.4 0.4 0.6],0.5,2); % all flys
            h.med{jj,mm}.Color(4) = 0.5;
            delete(h.patch)

            plot(HeadFree{amp(mm)}.CROSS.GrandMed.pat2head.delay(:,jj),...
                 HeadFree{amp(mm)}.CROSS.GrandMed.pat2head.maxCC(:,jj),'s','MarkerSize',10,...
                'MarkerEdgeColor','k','MarkerFaceColor',CList{amp(mm)})
            
            xlim([-1 1])
%             ylim([-20 20])
            if jj==HeadFree{amp(mm)}.n.Freq
                xlabel('Time Lag (s)','Interpreter','latex','FontSize',15)
            end
            if jj~=HeadFree{amp(mm)}.n.Freq
                xticks(0)
                xticklabels('')
            end
            ylabel('Head CC)','Interpreter','latex','FontSize',15)
    end
end
subplot(HeadFree{amp(1)}.n.Freq,1,1)
h.leg = legend([h.med{1,:}],legList);
title(h.leg,['Amplitude (' char(176) ')'])

% saveas(FIG,[root 'FIGURE\' filename '.fig']); % save .fig file
% print(gcf,[root 'FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
% disp('Saved to')
% disp(root)
end