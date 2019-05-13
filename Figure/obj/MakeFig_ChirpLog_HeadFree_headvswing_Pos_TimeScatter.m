function [FIG] = MakeFig_ChirpLog_HeadFree_headvswing_Pos_TimeScatter()
%% MakeFig_ChirpLog_HeadFree_headvswing_Pos_TimeScatter:
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

HeadFree = load(fullfile(root,FILE{1}),'TRIAL','U','N'); % load data structure

figNum = 1;
catIdx = [2 3];
xIdx = 1;

filename = 'ChirpLog_HeadFree_headvswing_Pos_TimeScatter'; % name of figure to save

HEAD = cell(HeadFree.N{1,3},1);
WING = HEAD;

% Store data by amplitude
pp = 1;
for jj = 1:HeadFree.N{1,3}
    for kk = 1:HeadFree.N{1,1}
        for ii = 1:size(HeadFree.TRIAL{kk,jj},1)
            HEAD{jj}(:,end+1) = HeadFree.TRIAL{kk,jj}{ii,catIdx(1)}.X(:,xIdx);
            WING{jj}(:,end+1) = HeadFree.TRIAL{kk,jj}{ii,catIdx(2)}.X(:,xIdx);
            pp = pp + 1; 
        end
    end
end

FIG = figure (figNum); clf % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1400 500];
FIG.Name = filename;
movegui(FIG,'center')
hold on

for jj = 1:HeadFree.N{1,3}
    ax = subplot(1,HeadFree.N{1,3},jj) ; hold on
    title([num2str(HeadFree.U{1,3}{1}(jj)) '$^{\circ}$'],'Interpreter','latex','FontSize',15)
    
    xData = HEAD{jj}(:);
    yData = WING{jj}(:);
    
  	[h] = scatplot(xData,yData);
    axis([-20 20 -8 8])
	delete(h.cb)
    
    xlabel('Head ($^{\circ}$)','Interpreter','latex','FontSize',15)
    if jj==1
        ylabel('Wing (V)','Interpreter','latex','FontSize',15)
    else
        yticks(0) ; yticklabels('')
    end

    % Calculate linear best fit
    [r,m,b] = regression(xData,yData,'one');
    text(20,7,['r =' num2str(r)]);
    xFit = linspace(-20,20,4000);
    yFit = m*xFit + b;
    plot(xFit,yFit,'r','LineWidth',5)   
end

saveas(FIG,[root 'FIGURE\' filename '.fig']); % save .fig file
% print(gcf,[root 'FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
end