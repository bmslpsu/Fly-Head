function [FIG] = MakeFig_ChirpLog_HeadFree_pat2head_Pos_TimeDelay()
%% MakeFig_ChirpLog_HeadFree_pat2head_Pos_TimeDelay:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'F:\DATA\Rigid_Data\';

% Select files
[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

HeadFree = load(fullfile(root,FILE{1}),'TRIAL','FLY','GRAND','U','N'); % load data structure

figNum = 1;
catIdx = 5;
xIdx = 1;

filename = 'ChirpLog_HeadFree_pat2head_Pos_TimeDelay'; % name of figure to save

% Get time difference
TD = cell(HeadFree.N{1,3},1);
for jj = 1:HeadFree.N{1,3}
    pp = 1;
    for kk = 1:HeadFree.N{1,1}
        for ii = 1:size(HeadFree.TRIAL{kk,jj},1)
            TD{jj}(pp,2) = HeadFree.TRIAL{kk,jj}{ii,catIdx}.TimeDiff(:,xIdx);
            TD{jj}(pp,1) = jj;
            pp = pp + 1;
        end
    end
end

TD = cat(1,TD{:});

FIG = figure (figNum); clf % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1000 400];
FIG.Name = filename;
movegui(FIG,'center')
hold on

h = boxplot(TD(:,2),TD(:,1),'Labels',{num2str(HeadFree.U{1,3}{1})},'Width',0.5,'Symbol','','Whisker',2);
ylim(0.1*[-1 1])
ylabel('Time Difference (s)','Interpreter','latex','FontSize',15)
xlabel('Amplitude ($^{\circ}$)','Interpreter','latex','FontSize',15)
% axis tight
ax = gca;
% set(gca,'FontSize',28)
hh = get(h(5,:),{'XData','YData'});
for kk = 1:size(hh,1)
   patch(hh{kk,1},hh{kk,2},'b');
end
set(findobj(FIG,'tag','Median'), 'Color', 'r');
set(findobj(FIG,'tag','Median'), 'LineWidth', 3);
set(findobj(FIG,'tag','Box'), 'Color', 'k');
set(findobj(FIG,'tag','Box'), 'LineWidth', 1);
set(findobj(FIG,'tag','Upper Whisker'), 'Color', 'k');
set(findobj(FIG,'tag','Upper Whisker'), 'LineWidth', 2);
set(findobj(FIG,'tag','Lower Whisker'), 'Color', 'k');
set(findobj(FIG,'tag','Lower Whisker'), 'LineWidth', 2);
ax.Children = ax.Children([end 1:end-1]);

saveas(FIG,[root 'FIGURE\' filename '.fig']); % save .fig file
% print(gcf,[root 'FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
end