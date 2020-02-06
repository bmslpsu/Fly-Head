function [] = MakeFig_Sine_HeadFree_pat_wing_Time()
%% MakeFig_Sine_HeadFree_pat_wing_Time:
%   INPUTS:
%       -
%   OUTPUTS:
%       -
%

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
FIG.Units = 'inches';
FIG.Position = [1 1 4 5];
movegui(FIG,'center')
FIG.Name = filename;
for ww = 1:nAmp
   FIG.Name = [FIG.Name '_' num2str(Amp(ww))];  
end

for ww = 1
    pp = 1;
    for jj = 1:HeadFree{ww}.N{1,3}
        freq = HeadFree{ww}.U{1,3}{1}(jj);
        
      	ax.R = subplot(HeadFree{ww}.N{1,3},1,pp) ; hold on
        h.stim = plot(HeadFree{ww}.TRIAL{5,jj}{2,1}.Time,HeadFree{ww}.TRIAL{5,jj}{2,1}.X(:,xIdx),...
            '-k','LineWidth',0.5);
        
        ax.R.YAxisLocation = 'right';
        ax.R.Color = 'none';
        ax.R.YColor = 'k';
        ax.R.YLabel.String = ['(' char(176) ')'];
        ax.R.YLabel.FontSize = 8;
        ax.R.XLim = [0 5];
        ax.R.YLim = 5*[-1 1];
        ax.R.YTick = 5*[-1 0 1];
        if pp==HeadFree{ww}.N{1,3}
            ax.R.XLabel.String = 'Time (s)';
            ax.R.XLabel.FontSize = 8;
        end
        
     	ax.L = axes; hold on
     	title([num2str(freq) ' Hz'],'FontSize',8,'FontWeight','bold')
        ax.L.Color = 'none';
    	ax.L.YLim = 2*[-1 1];
        ax.L.YColor = [0 0 0];
        ax.L.YLabel.String = '\Delta WBA (V)';
        ax.L.YLabel.FontSize = 8;
        ax.L.FontSize = 8;
        ax.L.XLim = [0 5];
        ax.L.XTickLabels = '';
        ax.L.XColor = 'none';
        
        h.patch = PlotPatch(HeadFree{ww}.GRAND{jj,catIdx}.Mean{2}{6}(:,xIdx),HeadFree{ww}.GRAND{jj,catIdx}.STD{2}{6}(:,xIdx),...
            HeadFree{ww}.GRAND{jj,catIdx}.Mean{2}{5},1,HeadFree{ww}.N{1,1},'r',[0.4 0.4 0.6],0.5,1);
        h.patch.EdgeColor = 'none';
        
        ax.L.Position = ax.R.Position;
%         uistack(ax.L,'top')
        
        pp = pp + 1;
    end
end

end