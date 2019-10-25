function [FIG] = MakeFig_Static_Head_Saccade_Pos_Vel_DIR()
%% MakeFig_Static_Head_Saccade_Pos_Vel_DIR:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

[FILE.Plus,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select positive saccade data', root, 'MultiSelect','off');

[FILE.Minus,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select negative saccade data', root, 'MultiSelect','off');

Plus    = load(fullfile(root,FILE.Plus),'SACCADE','U','N');
Minus      = load(fullfile(root,FILE.Minus),'SACCADE','U','N');

% clearvars -except SACCADE INTERVAL SACD Stim U N I TRIAL FLY GRAND CC Vel
CC = jet(Plus.N{1,3});
Wave = Plus.U{1,3}{1};

%% Saccade Position %%
clms = 6;
rows = ceil(Plus.N{1,3}/clms);

FIG = figure (1) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 clms*(4/3) rows*(3/2)];
FIG.Name = 'Saccade Position';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(Plus.N{1,3},1);
clear ax h
for jj = 1:Plus.N{1,3}
    ax(jj) = subplot(rows,clms,jj) ; hold on
    ax(jj).Title.String = [num2str(Wave(jj)) char(176)];
    ax(jj).Title.Color = 'k';
    
    cent = max(Plus.SACCADE.cIdx{jj},[],'all');
    span = (cent-5):(cent+5);
    
%     h.trial = plot(1000*Plus.SACCADE.Head.Time{jj}, Plus.SACCADE.Head.Position{jj}, 'Color', [0.7*CC(jj,:), 0.2]);

    PlotPatch(Plus.SACCADE.HeadStats.Position(jj).Median(span), Plus.SACCADE.HeadStats.Position(jj).STD(span), ...
        1000*Plus.SACCADE.HeadStats.Time(jj).Median(span), 1, 1, CC(jj,:), [0.7 0.7 0.7], 0.4, 3);
    
  	cent = max(Minus.SACCADE.cIdx{jj},[],'all');
    span = (cent-5):(cent+5);
    
%   	h.trial = plot(1000*Minus.SACCADE.Head.Time{jj}, Minus.SACCADE.Head.Position{jj}, 'Color', [0.5 0.5 0.5, 0.2]);

 	PlotPatch(Minus.SACCADE.HeadStats.Position(jj).Median(span), Minus.SACCADE.HeadStats.Position(jj).STD(span), ...
        1000*Minus.SACCADE.HeadStats.Time(jj).Median(span), 1, 1, CC(jj,:), [0.7 0.7 0.7], 0.4, 3);
end
set(ax,'FontSize',8,'Color','w','YColor','k','XColor','k','XLim',1000*0.05*[-1 1],'YLim',20*[-1 1])
XLabelHC = get(ax, 'XLabel');
set([XLabelHC{:}], 'String', 'Time (ms)')
YLabelHC = get(ax([1,4]), 'YLabel');
set([YLabelHC{:}], 'String', ['Head Position (' char(176) ')'])
linkaxes(ax(1:3),'y')
linkaxes(ax(4:6),'y')
linkaxes(ax,'x')

%% Saccade Velocity %%
FIG = figure (2) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 clms*(4/3) rows*(3/2)];
FIG.Name = 'Saccade Velocity';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
clear h
ax = gobjects(Plus.N{1,3},1);
for jj = 1:Plus.N{1,3}
    ax(jj) = subplot(rows,clms,jj) ; hold on
    ax(jj).Title.String = [num2str(Wave(jj)) char(176)];
    ax(jj).Title.Color = 'k';
    
    cent = max(Plus.SACCADE.cIdx{jj},[],'all');
    span = (cent-5):(cent+5);
    
%     h.trial = plot(1000*Plus.SACCADE.Head.Time{jj},Plus.SACCADE.Head.Velocity{jj}, 'Color', [0.7*CC(jj,:) , 0.2]);
    
	PlotPatch(Plus.SACCADE.HeadStats.Velocity(jj).Median(span), Plus.SACCADE.HeadStats.Velocity(jj).STD(span), ...
        1000*Plus.SACCADE.HeadStats.Time(jj).Median(span), 1, 1, CC(jj,:), [0.7 0.7 0.7], 0.4, 3);
    
 	cent = max(Minus.SACCADE.cIdx{jj},[],'all');
    span = (cent-5):(cent+5);
    
% 	h.trial = plot(1000*Minus.SACCADE.Head.Time{jj}, -Minus.SACCADE.Head.Velocity{jj}, 'Color', [0.5 0.5 0.5, 0.2]);
    
	PlotPatch(Minus.SACCADE.HeadStats.Velocity(jj).Median(span), Minus.SACCADE.HeadStats.Velocity(jj).STD(span), ...
        1000*Minus.SACCADE.HeadStats.Time(jj).Median(span), 1, 1, CC(jj,:), [0.7 0.7 0.7], 0.4, 3);
    
%     plot([-100 100],-sign(Vel(jj))*300*[1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',1)

%     if sign(Vel(jj))==1
%         ax(jj).YLim = 1100*[-1 0.2];
%     else
%         ax(jj).YLim =  1100*[-0.2 1];
%   	end
    
    if sign(Wave(jj))==1
        ax(jj).YLim = 1100*[-1 0.2];
    else
        ax(jj).YLim =  1100*[-0.2 1];
    end
    
end
set(ax,'FontSize',8,'Color','w','YColor','k','XColor','k','XLim',1000*0.05*[-1 1])
XLabelHC = get(ax, 'XLabel');
set([XLabelHC{:}], 'String', 'Time (ms)')
YLabelHC = get(ax([1,4]), 'YLabel');
set([YLabelHC{:}], 'String', ['Head Velocity (' char(176) '/s)'])
linkaxes(ax(1:3),'y')
linkaxes(ax(4:6),'y')
linkaxes(ax,'x')

end