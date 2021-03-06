function [FIG] = MakeFig_Ramp_Head_Removed()
%% MakeFig_Ramp_Head_Removed:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';
[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

% load(fullfile(root,FILE{1}),'SACCADE','INTERVAL','SACD','Stim','U','N','I','TRIAL','FLY','GRAND');
load(fullfile(root,FILE{1}),'SACD','Stim','U','N','TRIAL','FLY','GRAND');

clearvars -except SACCADE INTERVAL SACD Stim U N I TRIAL FLY GRAND

CC = repmat({[1 0 0],[0 1 0],[0 0 1]},1,2);
Vel = 3.75*U{1,3}{1};

%% Saccade Removed Time Domain %%
FIG = figure (1) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 6 4];
FIG.Name = 'Removed Saccades';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
hold on
ax = gobjects(N{1,3},1);
rmvIdx = 4;
for jj = 1:N.vel
    ax(jj) = subplot(2,3,jj) ; hold on
    ax(jj).Title.String = [num2str(Vel(jj)) ' (' char(176) '/s)'];
    ax(jj).XLabel.String = 'Time(s)';
    ax(jj).YLabel.String = ['Head Angle (' char(176) ')'];
    for kk = 1:N.fly
        for ii = 1:size(TRIAL{kk,jj},1)
            h.trial = plot(TRIAL{kk,jj}{ii,rmvIdx}.Time, TRIAL{kk,jj}{ii,rmvIdx}.X(:,1), 'Color', 0.5*CC{jj});
            h.trial.Color(4) = 0.2;
        end
        plot(FLY{jj}{kk,rmvIdx}.Mean{5}, FLY{jj}{kk,rmvIdx}.Mean{6}(:,1), 'Color', 0.7*CC{jj},'LineWidth',1)
    end
    plot(GRAND{jj,rmvIdx}.Mean{1}{5}, GRAND{jj,rmvIdx}.Mean{1}{6}(:,1), 'Color', CC{jj},'LineWidth',3)
    plot(GRAND{jj,rmvIdx}.Mean{1}{5},Stim(:,jj), '--', 'Color', 'k', 'LineWidth',2)
end
set(ax,'XLim',[0 10])
set(ax(1:3),'YLim',200*[-0.1 1])
set(ax(4:6),'YLim',200*[-1 0.1])
set(ax([1:4,6]),'XColor','none')
set(ax([2:3,5:6]),'YColor','none')
set(ax,'FontSize',8)

%% Saccade Removed Time Domain ALL %%
FIG = figure (2) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 6 4];
FIG.Name = 'Removed Saccades';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
hold on
clear ax h
ax(1) = subplot(1,1,1) ; hold on
ax(1).XLabel.String = 'Time(s)';
ax(1).YLabel.String = ['Head Angle (' char(176) ')'];
rmvIdx = 4;
xIdx = 1;
for jj = 1:N.vel
%     for kk = 1:N.fly
%         for ii = 1:size(TRIAL{kk,jj},1)
%             h.trial = plot(TRIAL{kk,jj}{ii,rmvIdx}.Time, TRIAL{kk,jj}{ii,rmvIdx}.X(:,1), 'Color', 0.5*CC{jj});
%             h.trial.Color(4) = 0.2;
%         end
%         plot(FLY{jj}{kk,rmvIdx}.Mean{5}, FLY{jj}{kk,rmvIdx}.Mean{6}(:,1), 'Color', 0.7*CC{jj},'LineWidth',1)
%     end
    % h.mean(jj) = plot(GRAND{jj,rmvIdx}.Mean{1}{5}, GRAND{jj,rmvIdx}.Mean{1}{6}(:,1), 'Color', CC{jj},'LineWidth',3);
    
	[h.std(jj),h.mean(jj)] = PlotPatch(GRAND{jj,rmvIdx}.Mean{1}{6}(:,xIdx), GRAND{jj,rmvIdx}.STD{1}{6}(:,xIdx), ...
        GRAND{jj,rmvIdx}.Mean{1}{5}, 1, 1, CC{jj}, 0.5*CC{jj}, 0.2, 3);
    
    
%     plot(GRAND{jj,rmvIdx}.Mean{1}{5},Stim(:,jj), '--', 'Color', 'k', 'LineWidth',2)
end
set(ax,'XLim',[0 10])
% set(ax,'YLim',200*[-1 1])
% set(ax,'XColor','none')
% set(ax,'YColor','none')
set(ax,'FontSize',8)
uistack(h.std,'top')
uistack(h.mean,'top')

%% Saccade Removed Frequency Domain %%
FIG = figure (2) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 6 4];
FIG.Name = 'Removed Saccades';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
hold on
ax = gobjects(N{1,3},1);
rmvIdx = 2;
for jj = 1:N.vel
    ax(jj) = subplot(2,3,jj) ; hold on
    ax(jj).Title.String = [num2str(Vel(jj)) ' (' char(176) '/s)'];
    ax(jj).XLabel.String = 'Frequency (Hz)';
    ax(jj).YLabel.String = ['Head Angle (' char(176) ')'];
    for kk = 1:N.fly
        for ii = 1:size(TRIAL{kk,jj},1)
            h.trial = plot(TRIAL{kk,jj}{ii,rmvIdx}.Fv, TRIAL{kk,jj}{ii,rmvIdx}.Mag(:,1), 'Color', 0.5*CC{jj});
            h.trial.Color(4) = 0.2;
        end
        plot(FLY{jj}{kk,rmvIdx}.Mean{7}, FLY{jj}{kk,rmvIdx}.Mean{8}(:,1), 'Color', 0.7*CC{jj},'LineWidth',1)
    end
    plot(GRAND{jj,rmvIdx}.Mean{1}{7}, GRAND{jj,rmvIdx}.Mean{1}{8}(:,1), 'Color', CC{jj},'LineWidth',3)
end
set(ax,'XLim',[0 40])
set(ax,'YLim',[0 1])
% set(ax([1:4,6]),'XColor','none')
set(ax([2:3,5:6]),'YColor','none')
set(ax,'FontSize',8)

%% Saccade Removed Frequency Domain ALL %%
FIG = figure (2) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 3 2];
FIG.Name = 'Removed Saccades';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
hold on
clear ax
ax = subplot(1,1,1) ; hold on
rmvIdx = 2;
for jj = 1:N.vel
%     ax.Title.String = [num2str(Vel(jj)) ' (' char(176) '/s)'];
    ax.XLabel.String = 'Frequency (Hz)';
    ax.YLabel.String = ['Head Angle (' char(176) ')'];
%     for kk = 1:N.fly
%         for ii = 1:size(TRIAL{kk,jj},1)
%             h.trial = plot(TRIAL{kk,jj}{ii,rmvIdx}.Fv, TRIAL{kk,jj}{ii,rmvIdx}.Mag(:,1), 'Color', 0.5*CC{jj});
%             h.trial.Color(4) = 0.3;
%         end
% %         plot(FLY{jj}{kk,rmvIdx}.Mean{7}, FLY{jj}{kk,rmvIdx}.Mean{8}(:,1), 'Color', 0.7*CC{jj},'LineWidth',1)
%     end
    plot(GRAND{jj,rmvIdx}.Mean{1}{7}, GRAND{jj,rmvIdx}.Mean{1}{8}(:,1), 'Color', CC{jj},'LineWidth',3)
    
%     PlotPatch(GRAND{jj,rmvIdx}.Mean{1}{8}(:,1), GRAND{jj,rmvIdx}.STD{1}{8}(:,1), ...
%         GRAND{jj,rmvIdx}.Mean{1}{7}, 2, N.fly, CC{jj}, [0.7 0.7 0.7], 0.4, 3);
end
set(ax,'XLim',[5 30])
set(ax,'YLim',[0 0.7])
% set(ax([1:4,6]),'XColor','none')
set(ax,'FontSize',8)

end