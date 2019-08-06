function [] = MakeFig_Static_Correlation()
%% MakeFig_Static_Correlation:
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

clearvars -except HeadFree Amp nAmp

filename = 'Static_Correlation';

headIdx = 2;
wingIdx = 3;
xIdx = 1;

nFreq = HeadFree{1}.N{1,3};

HEAD = cell(nFreq,nAmp);
WING = cell(nFreq,nAmp);
for ww = 1:nAmp
    nFly = HeadFree{ww}.N{1,1};
    for jj = 1:nFreq
        pp = 1;
        for kk = 1:nFly
            for ii = 1:size(HeadFree{ww}.TRIAL{kk,jj},1) % trials
                HEAD{jj,ww}(:,pp) = HeadFree{ww}.TRIAL{kk,jj}{ii,headIdx}.X(:,xIdx);
                WING{jj,ww}(:,pp) = HeadFree{ww}.TRIAL{kk,jj}{ii,wingIdx}.X(:,xIdx);
                pp =  pp + 1;
            end
        end
    end
end
%%

FIG = figure (1); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 3 3];
FIG.Name = filename;
movegui(FIG,'center')

for ww = 1:nAmp
    pp = 1;
    for jj = 1
        ax = subplot(1,1,pp) ; hold on
        ax.Title.FontWeight = 'bold';
      	ax.FontSize = 8;
        ax.YLabel.FontSize = 8;
        ax.XLim = 20*[-1 1];
        ax.YLim = 6*[-1 1];
        ax.XLabel.FontSize = 8;
        ax.YLabel.FontSize = 8;
        
        ax.XLabel.String = ['Head (' char(176) ')'];
        ax.YLabel.String = 'Wing (V)';
        
       	xData = HEAD{jj,ww};
        yData = WING{jj,ww};
        
     	[h] = scatplot(xData(:),yData(:));
        delete(h.cb)
        h.fit = lsline;
        
        % Calculate linear best fit
        [r,m,b] = regression(xData',yData','one');
        text(mean(ax.XLim),mean(ax.YLim),['r =' num2str(r)]);
        xFit = linspace(-20,20,4000);
        yFit = m*xFit + b;
        plot(xFit,yFit,'r','LineWidth',2)
        
        pp =  pp + 1;
    end
end
%%
FIG = figure (2); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 3 3];
FIG.Name = filename;
movegui(FIG,'center')
for ww = 1:nAmp
    pp = 1;
    for jj = 1
        ax = subplot(1,1,pp) ; hold on
        ax.Title.FontWeight = 'bold';
      	ax.FontSize = 8;
        ax.YLabel.FontSize = 8;
        ax.XLim = 20*[-1 1];
        ax.YLim = 6*[-1 1];
        ax.XLabel.FontSize = 8;
        ax.YLabel.FontSize = 8;
        
        ax.XLabel.String = ['Head (' char(176) ')'];
        ax.YLabel.String = 'Wing (V)';
        
       	xData = HEAD{jj,ww};
        yData = WING{jj,ww};
        
        cmap = jet(size(xData,2));

        for kk = 1:size(xData,2)
            h.scat = scatter(xData(:,kk),yData(:,kk),2,'MarkerFaceColor',cmap(kk,:),'MarkerFaceAlpha',0.2,...
                                'MarkerEdgeColor','none');
        end
        
        h.fit = lsline;
        
        % Calculate linear best fit
        [r,m,b] = regression(xData,yData,'one');
        text(mean(ax.XLim),mean(ax.YLim),['r =' num2str(r)]);
        xFit = linspace(-20,20,4000);
        yFit = m*xFit + b;
        plot(xFit,yFit,'r','LineWidth',2)
        
        pp =  pp + 1;
    end
end
end