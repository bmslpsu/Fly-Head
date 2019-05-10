function [] = MakeFig_Static_HeadFree_Head_Freq()
%% MakeFig_Ramp_HeadFree_Head_Freq: head frequency spectrum plot for ramp
%   INPUTS:
%       -
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE INPUT %
root = 'F:\DATA\Rigid_Data\';
figNum = 1;
%---------------------------------------------------------------------------------------------------------------------------------
filename = 'Static_HeadFree_Head_Freq'; % name of figure to save
HeadFree = load(fullfile(root, 'Static_HeadFree_DATA.mat'),'TRIAL','FLY','GRAND','U','N'); % load data structure

FIG = figure (figNum); % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1100 800];
FIG.Name = filename;
hold on
catIdx = 1;
xIdx = 1;

% Trials
pp = 1;
for jj = 1:HeadFree.N{1,3} % spatial frequencies
    subplot(HeadFree.N{1,3},1,pp) ; hold on ; xlim([0 30]) ; ylim([0 1])
    for kk = 1:HeadFree.N{1,1} % flys
        for ii = 1: size(HeadFree.TRIAL{kk,jj},1) % trials
            h.Trial = plot(HeadFree.TRIAL{kk,jj}{ii,catIdx}.Fv(:,xIdx),HeadFree.TRIAL{kk,jj}{ii,catIdx}.Mag(:,xIdx),...
                '-','Color',[0.5 0.5 0.5 0.25],'LineWidth',1);
        end
    end
    pp = pp + 1;
end


% Fly Stats
pp = 1;
for jj = 1:HeadFree.N{1,3} % spatial frequencies
    subplot(HeadFree.N{1,3},1,pp) ; hold on ; xlim([0 30]) ; ylim([0 1])
    for kk = 1:HeadFree.N{1,1} % flys
        h.Fly = plot(HeadFree.FLY{jj}{kk,catIdx}.Mean{7}(:,xIdx),HeadFree.FLY{jj}{kk,catIdx}.Mean{8}(:,xIdx),'-','LineWidth',2);
        h.Fly.Color(4) = 0.5;
    end
    pp = pp + 1;
end


% Grand Stats
vel = HeadFree.U{1,4}{1};
pp = 1;
for jj = 1:HeadFree.N{1,3} % spatial frequencies
    subplot(HeadFree.N{1,3},1,pp) ; hold on ; xlim([0 20]) ; ylim([0 1])
%         grid on ; grid minor
    freq = HeadFree.U{1,3}{1}(jj);
    
    if pp==1
        title([num2str(vel) char(176) '/s'],'FontSize',15,'FontWeight','bold')
    end
    
   	ylabel([num2str(freq) char(176)],'FontSize',15,'FontWeight','bold')
    
    if pp==HeadFree.N{1,3}
        xlabel('Frequency (Hz)','FontSize',12,'FontWeight','normal')
    else
        xticks('')
    end

    h.patch = PlotPatch(HeadFree.GRAND{jj,catIdx}.Mean{2}{8}(:,xIdx),HeadFree.GRAND{jj,catIdx}.STD{2}{8}(:,xIdx),...
        HeadFree.GRAND{jj,catIdx}.Mean{2}{7}(:,xIdx),2,HeadFree.N{1,1},'k',[0.4 0.4 0.6],0.5,3);

    pp = pp + 1;
end

saveas(FIG,['F:\DATA\Rigid_Data\FIGURE\' filename '.fig']); % save .fig file
print(FIG,['F:\DATA\Rigid_Data\FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
end