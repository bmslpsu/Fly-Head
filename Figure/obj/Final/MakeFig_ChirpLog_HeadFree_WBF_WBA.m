function [FIG] = MakeFig_ChirpLog_HeadFree_WBF_WBA()
%% MakeFig_ChirpLog_HeadFree_WBF_WBA:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

[FILE,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
FILE = cellstr(FILE)';

HeadFree = load(fullfile(root,FILE{1}),'TRIAL','GRAND','U','N');

clearvars -except HeadFree Amp nAmp

catIdx = 3;

filename = 'MakeFig_ChirpLog_HeadFree_WBF_WBA';

%% WBF
FIG = figure (1); clf
FIG.Color = 'w';
FIG.Position = [100 100 800 400];
FIG.Name = filename;
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')

ax = gca; hold on
ax.FontSize = 12;
ax.XLabel.String = 'Time (s)';
ax.XLabel.FontSize = 14;
ax.YLabel.String = 'WBF (Hz)';
ax.YLabel.FontSize = 14;
ax.YLim = [200 280];

cList = prism(HeadFree.N.Amp);
clear h
for jj = 1:HeadFree.N.Amp
    for kk = 1:HeadFree.N.fly
        for ii = 1:size(HeadFree.TRIAL{kk,jj},1)
            h.trial = plot(HeadFree.TRIAL{kk,jj}{ii,catIdx}.Time, HeadFree.TRIAL{kk,jj}{ii,catIdx}.WBF,...
                'Color',cList(jj,:),'LineWidth',0.5);
            h.trial.Color(4) = 0.2;
        end
    end
	h.grand(jj) = plot(HeadFree.GRAND{jj,catIdx}.Mean{1}{5}, HeadFree.GRAND{jj,catIdx}.Mean{1}{26},...
        'Color','k','LineWidth',3);
end
uistack(h.grand,'top')

%% WBA
FIG = figure (2); clf
FIG.Color = 'w';
FIG.Position = [100 100 800 400];
FIG.Name = filename;
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')

ax = gca; hold on
ax.FontSize = 12;
ax.XLabel.String = 'Time (s)';
ax.XLabel.FontSize = 14;
ax.YLabel.String = 'WBA (Hz)';
ax.YLabel.FontSize = 14;
% ax.YLim = [200 280];

cList = prism(HeadFree.N.Amp);
clear h
for jj = 1:HeadFree.N.Amp
    for kk = 1:HeadFree.N.fly
        for ii = 1:size(HeadFree.TRIAL{kk,jj},1)
            h.trial = plot(HeadFree.TRIAL{kk,jj}{ii,catIdx}.Time, HeadFree.TRIAL{kk,jj}{ii,catIdx}.WBA(:,3),...
                'Color',cList(jj,:),'LineWidth',0.5);
            h.trial.Color(4) = 0.2;
        end
    end
	h.grand(jj) = plot(HeadFree.GRAND{jj,catIdx}.Mean{1}{5}, HeadFree.GRAND{jj,catIdx}.Mean{1}{27}(:,3),...
        'Color','k','LineWidth',3);
end
uistack(h.grand,'top')
end