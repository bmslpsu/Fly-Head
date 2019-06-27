function [] = MakeFig_Sine_HeadFree_pat_wing_Time()
%% MakeFig_Sine_HeadFree_pat_wing_Time:
%   INPUTS:
%       -
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

[FILES,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','on');
FILES = cellstr(FILES)';

nAmp = length(FILES);
Amp = nan(nAmp,1);
for ww = 1:nAmp
    filedata = textscan(FILES{ww}, '%s', 'delimiter', '_');
    Amp(ww) = str2double(filedata{1}{3});
end

HeadFree = cell(nAmp,1);
for ww = 1:nAmp
    HeadFree{ww} = load(fullfile(root,FILES{ww}),'TRIAL','GRAND','U','N');
end
%%
clearvars -except HeadFree Amp nAmp

filename = ['Sine_HeadFree_pat_wing_Time_' num2str(1)];

hold on
catIdx = 3;
xIdx = 1;
figNum = 1;

FIG = figure (figNum); clf
FIG.Color = 'w';
FIG.Position = [100 100 600 700];
movegui(FIG,'center')
FIG.Name = filename;
for ww = 1:nAmp
   FIG.Name = [FIG.Name '_' num2str(Amp(ww))];  
end

for ww = 4
    pp = 1;
    for jj = 1:HeadFree{ww}.N{1,3}
        freq = HeadFree{ww}.U{1,3}{1}(jj);
     	subplot(HeadFree{ww}.N{1,3},1,pp) ; hold on ; xlim([0 10]) ; 
     	title([num2str(freq) ' Hz'],'FontSize',14,'FontWeight','bold')
        
        yyaxis left ; ax.L = gca;
        ax.L.YLim = 2*[-1 1];
        ax.L.YColor = [0 0 0];
        ax.L.YLabel.String = ['Wing (V)'];
        ax.L.YLabel.FontSize = 14;
        ax.L.FontSize = 12;
        ax.L.XLim = [0 5];
        if pp==HeadFree{ww}.N{1,3}
            ax.L.XLabel.String = 'Time (s)';
            ax.L.XLabel.FontSize = 14;
        else
            ax.L.XTickLabels = [];
        end
        
        h.patch = PlotPatch(HeadFree{ww}.GRAND{jj,catIdx}.Mean{2}{6}(:,xIdx),HeadFree{ww}.GRAND{jj,catIdx}.STD{2}{6}(:,xIdx),...
            HeadFree{ww}.GRAND{jj,catIdx}.Mean{2}{5},2,HeadFree{ww}.N{1,1},'r',[0.4 0.4 0.6],0.5,2);
        
        yyaxis right ; ax.R = gca;
        ax.R.YColor = [0 1 0];
        ax.R.YLabel.String = ['Stimulus (' char(176) ')'];
        ax.R.YLabel.FontSize = 14;
        ax.R.XLim = [0 5];
        if pp==HeadFree{ww}.N{1,3}
            ax.R.XLabel.String = 'Time (s)';
            ax.R.XLabel.FontSize = 14;
        end
        plot(HeadFree{ww}.TRIAL{5,jj}{2,1}.Time,HeadFree{ww}.TRIAL{5,jj}{2,1}.X(:,xIdx),'-g','LineWidth',1.5)
        
        pp = pp + 1;
    end
end

end