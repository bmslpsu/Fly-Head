function [FIG] = MakeFig_Ramp_Head_Saccade_Pos_Vel_Anti_CO()
%% MakeFig_Ramp_Head_Saccade_Pos_Vel_Anti_CO:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

[FILE.Anti,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select anti-directional data', root, 'MultiSelect','off');

[FILE.Co,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select co-directional data', root, 'MultiSelect','off');

Anti    = load(fullfile(root,FILE.Anti),'SACCADE','U','N');
Co      = load(fullfile(root,FILE.Co),'SACCADE','U','N');

CC = repmat({[1 0 0],[0 1 0],[0 0 1]},1,2);

Vel = 3.75*Anti.U{1,3}{1};
% clearvars -except SACCADE INTERVAL SACD Stim U N I TRIAL FLY GRAND CC Vel

%% Saccade Position %%
FIG = figure (1) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 4 3];
FIG.Name = 'Saccade Position';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(Anti.N{1,3},1);
clear ax h
for jj = 1:Anti.N{1,3}
    ax(jj) = subplot(ceil(Anti.N{1,3}/3),3,jj) ; hold on
    ax(jj).Title.String = [num2str(Vel(jj)) ' (' char(176) '/s)'];
    ax(jj).Title.Color = 'k';
    
    cent = max(Anti.SACCADE.cIdx{jj},[],'all');
    span = (cent-5):(cent+5);
    
%     h.trial = plot(1000*Anti.SACCADE.Head.Time{jj}, Anti.SACCADE.Head.Position{jj}, 'Color', [0.7*CC{jj}, 0.2]);

    PlotPatch(Anti.SACCADE.HeadStats.Position(jj).Median(span), Anti.SACCADE.HeadStats.Position(jj).STD(span), ...
        1000*Anti.SACCADE.HeadStats.Time(jj).Median(span), 1, 1, CC{jj}, [0.7 0.7 0.7], 0.4, 3);
    
  	cent = max(Co.SACCADE.cIdx{jj},[],'all');
    span = (cent-5):(cent+5);
    
%   	h.trial = plot(1000*Co.SACCADE.Head.Time{jj}, Co.SACCADE.Head.Position{jj}, 'Color', [0.5 0.5 0.5, 0.2]);

 	PlotPatch(-Co.SACCADE.HeadStats.Position(jj).Median(span), Co.SACCADE.HeadStats.Position(jj).STD(span), ...
        1000*Co.SACCADE.HeadStats.Time(jj).Median(span), 1, 1, [0.5 0.5 0.5], [0.7 0.7 0.7], 0.4, 3);
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
ax = gobjects(Anti.N{1,3},1);
for jj = 1:Anti.N{1,3}
    ax(jj) = subplot(ceil(Anti.N{1,3}/3),3,jj) ; hold on
    ax(jj).Title.String = [num2str(Vel(jj)) ' (' char(176) '/s)'];
    ax(jj).Title.Color = 'k';
    
    cent = max(Anti.SACCADE.cIdx{jj},[],'all');
    span = (cent-5):(cent+5);
    
%     h.trial = plot(1000*Anti.SACCADE.Head.Time{jj},Anti.SACCADE.Head.Velocity{jj}, 'Color', [0.7*CC{jj} , 0.2]);
    
	PlotPatch(Anti.SACCADE.HeadStats.Velocity(jj).Median(span), Anti.SACCADE.HeadStats.Velocity(jj).STD(span), ...
        1000*Anti.SACCADE.HeadStats.Time(jj).Median(span), 1, 1, CC{jj}, [0.7 0.7 0.7], 0.4, 3);
    
 	cent = max(Co.SACCADE.cIdx{jj},[],'all');
    span = (cent-5):(cent+5);
    
% 	h.trial = plot(1000*Co.SACCADE.Head.Time{jj}, -Co.SACCADE.Head.Velocity{jj}, 'Color', [0.5 0.5 0.5, 0.2]);
    
	PlotPatch(-Co.SACCADE.HeadStats.Velocity(jj).Median(span), Co.SACCADE.HeadStats.Velocity(jj).STD(span), ...
        1000*Co.SACCADE.HeadStats.Time(jj).Median(span), 1, 1, [0.5 0.5 0.5], [0.7 0.7 0.7], 0.4, 3);
    
%     plot([-100 100],-sign(Vel(jj))*300*[1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',1)

%     if sign(Vel(jj))==1
%         ax(jj).YLim = 1100*[-1 0.2];
%     else
%         ax(jj).YLim =  1100*[-0.2 1];
%   	end
    
    if sign(Vel(jj))==1
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