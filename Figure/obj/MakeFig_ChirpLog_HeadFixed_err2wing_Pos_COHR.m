function [FIG] = MakeFig_ChirpLog_HeadFixed_err2wing_Pos_COHR()
%% MakeFig_ChirpLog_HeadFixed_err2wing_Pos_COHR:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

% Select files
[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

HeadFree = load(fullfile(root,FILE{1}),'TRIAL','FLY','GRAND','U','N'); % load data structure

figNum = 1;
catIdx = 3;
xIdx = 1;

filename = 'ChirpLog_HeadFixed_err2wing_Pos_COHR'; % name of figure to save

FIG = figure (figNum); clf % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1100 400];
FIG.Name = filename;
movegui(FIG,'center')
hold on

% Trials
pp = 1;
for jj = 1:HeadFree.N{1,3}
    for kk = 1:HeadFree.N{1,1}
        for ii = 1:size(HeadFree.TRIAL{kk,jj},1)
            subplot(1,HeadFree.N{1,3},pp) ; hold on ; xlim([0 12.]) ; ylim([0 0.5])
                h.Trial = plot(HeadFree.TRIAL{kk,jj}{ii,catIdx}.CoherenceFV(:,xIdx),HeadFree.TRIAL{kk,jj}{ii,catIdx}.Coherence(:,xIdx),...
                    '-','Color',[0.5 0.5 0.5 0.1],'LineWidth',1);
        end
    end
    pp = pp + 1;
end

% Fly Stats
pp = 1;
for jj = 1:HeadFree.N{1,3}
    for kk = 1:HeadFree.N{1,1}
        subplot(1,HeadFree.N{1,3},pp) ; hold on ; xlim([0 12]) ; ylim([0 1])
            h.Fly = plot(HeadFree.FLY{jj}{kk,catIdx}.Mean{8}(:,xIdx),HeadFree.FLY{jj}{kk,catIdx}.Median{7}(:,xIdx),...
                '-','LineWidth',1);
            h.Fly.Color(4) = 0.5;
    end
    pp = pp + 1;
end

% Grand Stats
pp = 1;
for jj = 1:HeadFree.N{1,3}
    ax = subplot(1,HeadFree.N{1,3},pp) ; hold on ; xlim([0.1 12]) ; ylim([0 1])
    title([num2str(HeadFree.U{1,3}{1}(jj)) char(176)],'FontSize',15)
        h.patch = PlotPatch(HeadFree.GRAND{jj,catIdx}.Mean{2}{7}(:,xIdx),HeadFree.GRAND{jj,catIdx}.STD{2}{7}(:,xIdx),...
            HeadFree.GRAND{jj,catIdx}.Mean{2}{8}(:,xIdx),3,HeadFree.N{1,1},'k',[0.4 0.4 0.6],0.5,3);
               
        if pp==1
            ylabel('Error to Wing Coherence','Interpreter','latex','FontSize',15)
       	else
            yticks('')
        end
        
        xlabel('Frequency (Hz)','Interpreter','latex','FontSize',15)
        
        ax.XTick = sort([0.1 ax.XTick]);

    pp = pp + 1;
end

saveas(FIG,[root 'FIGURE\' filename '.fig']); % save .fig file
% print(gcf,[root 'FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
end