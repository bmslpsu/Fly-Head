function [] = MakeFig_Asymmetry_HeadFree_Head_Freq()
%% MakeFig_Asymmetry_HeadFree_Head_Freq: head frequency spectrum plot for ramp
%   INPUTS:
%       -
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE INPUT %
root = 'H:\EXPERIMENTS\Experiment_Asymmetry_Verification\DATA\';
figNum = 1;
%---------------------------------------------------------------------------------------------------------------------------------
filename = 'Asymmetry_HeadFree_Head_Freq'; % name of figure to save
HeadFree = load([root 'Asymmetry_HeadFree_DATA.mat'],...
    'TRIAL','FLY','GRAND','U','N'); % load data structure

close all
FIG = figure (figNum); % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1100 800];
FIG.Name = filename;
hold on
catIdx = 2;
xIdx = 1;

% Trials
hold on ; xlim([0 30]) ; ylim([0 1])
for kk = 1:HeadFree.N{1,1} % flys
    for ii = 1: size(HeadFree.TRIAL{kk,1},1) % trials
        h.Trial = plot(HeadFree.TRIAL{kk,1}{ii,catIdx}.Fv(:,xIdx),HeadFree.TRIAL{kk,1}{ii,catIdx}.Mag(:,xIdx),...
            '-','Color',[0.5 0.5 0.5 0.25],'LineWidth',1);
    end
end

% Fly Stats
hold on ; xlim([0 30]) ; ylim([0 1])
for kk = 1:HeadFree.N{1,1} % flys
    h.Fly = plot(HeadFree.FLY{1}{kk,catIdx}.Mean{7}(:,xIdx),HeadFree.FLY{1}{kk,catIdx}.Mean{8}(:,xIdx),'-','LineWidth',2);
    h.Fly.Color(4) = 0.5;
end

% Grand Stats
hold on ; xlim([0 40]) ; ylim([0 1]) ; grid on ; grid minor
vel = 45;
freq = 0;
tmpFreq = vel/freq;
title([num2str(vel) char(176) '/s'],'FontSize',15,'FontWeight','bold')
ylabel([num2str(freq) char(176)],'FontSize',15,'FontWeight','bold')
xlabel('Frequency (Hz)','FontSize',12,'FontWeight','normal')
text(15,0.9,[num2str(tmpFreq) ' cyc/s'])

h.patch = PlotPatch(HeadFree.GRAND{1,catIdx}.Mean{2}{8}(:,xIdx),HeadFree.GRAND{1,catIdx}.STD{2}{8}(:,xIdx),...
    HeadFree.GRAND{1,catIdx}.Mean{2}{7}(:,xIdx),2,HeadFree.N{1,1},'k',[0.4 0.4 0.6],0.5,3);

saveas(FIG,[root 'FIGURE\' filename '.fig']); % save .fig file
print(gcf,[root 'FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
disp('Saved to')
disp(root)
end