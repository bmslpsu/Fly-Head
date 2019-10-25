function [FIG] = MakeFig_Ramp_Head_Saccade_Pos_Vel()
%% MakeFig_Ramp_Head_Saccade_Pos_Vel:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

load(fullfile(root,FILE{1}),'SACCADE','U','N');

CC = repmat({[1 0 0],[0 1 0],[0 0 1]},1,2);

Vel = 3.75*U{1,3}{1};
clearvars -except SACCADE INTERVAL SACD Stim U N I TRIAL FLY GRAND CC Vel

%% Saccade Position %%
FIG = figure (1) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 4 3];
FIG.Name = 'Saccade Position';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(N{1,3},1);
clear ax h
for jj = 1:N{1,3}
    ax(jj) = subplot(ceil(N{1,3}/3),3,jj) ; hold on
    ax(jj).Title.String = [num2str(Vel(jj)) ' (' char(176) '/s)'];
    ax(jj).Title.Color = 'k';
    
    cent = max(SACCADE.cIdx{jj},[],'all');
    span = (cent-5):(cent+5);
    
    h.trial = plot(1000*SACCADE.Head.Time{jj}, SACCADE.Head.Position{jj}, 'Color', [0.7*CC{jj} , 0.2]);
    
    PlotPatch(SACCADE.HeadStats.Position(jj).Median(span), SACCADE.HeadStats.Position(jj).STD(span), ...
        1000*SACCADE.HeadStats.Time(jj).Median(span), 1, 1, CC{jj}, [0.7 0.7 0.7], 0.4, 3);
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
FIG.Position = [2 2 4 3];
FIG.Name = 'Saccade Velocity';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
clear h
ax = gobjects(N{1,3},1);
for jj = 1:N{1,3}
    ax(jj) = subplot(ceil(N{1,3}/3),3,jj) ; hold on
    ax(jj).Title.String = [num2str(Vel(jj)) ' (' char(176) '/s)'];
    ax(jj).Title.Color = 'k';
    
    cent = max(SACCADE.cIdx{jj},[],'all');
    span = (cent-5):(cent+5);
    
    h.trial = plot(1000*SACCADE.Head.Time{jj},SACCADE.Head.Velocity{jj}, 'Color', [0.7*CC{jj} , 0.2]);
    
	PlotPatch(SACCADE.HeadStats.Velocity(jj).Median(span), SACCADE.HeadStats.Velocity(jj).STD(span), ...
        1000*SACCADE.HeadStats.Time(jj).Median(span), 1, 1, CC{jj}, [0.7 0.7 0.7], 0.4, 3);
    
%     plot([-100 100],-sign(Vel(jj))*300*[1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',1)
    
%     if sign(Vel(jj))==1
%         ax(jj).YLim = 1100*[-1 0.2];
%     else
%         ax(jj).YLim =  1100*[-0.2 1];
%   	end
    
    if sign(Vel(jj))==1
        ax(jj).YLim = 1100*[-0.2 1];
    else
        ax(jj).YLim =  1100*[-1 0.2];
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