function [] = MakeFig_Static_Correlation()
%% MakeFig_Static_Correlation:
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

% clearvars -except HeadFree Amp nAmp

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
FIG.Position = [2 2 3*3 2*3];
FIG.Name = filename;
movegui(FIG,'center')

ax = gobjects(nAmp,nFreq);
pp = 1;
for ww = 1:nAmp
    for jj = [2:4,[6 5 1]]
       	xData = HEAD{jj,ww};
        yData = WING{jj,ww};
        
        ax(pp) = subplot(2,3,pp) ; hold on
            ax(pp).Title.String = num2str(HeadFree{ww}.U{1,3}{1}(jj));
            xlabel('Head (°)')
            ylabel('Wing (V)')
        
            [h] = scatplot(xData(:),yData(:));
            delete(h.cb)
            % h.fit = lsline;
        
            % Calculate linear best fit
            [r,m,b] = regression(xData',yData','one');
            text(ax(pp).XLim(1) + 12, ax(pp).YLim(1) + 0.3,['r =' num2str(r)]);
            xFit = linspace(-20,20,4000);
            yFit = m*xFit + b;
            plot(xFit,yFit,'r','LineWidth',2)
        
        pp =  pp + 1;
    end
end

set(ax, 'FontSize', 8, 'LineWidth', 1.5, 'XLim', 20*[-1 1], 'YLim', 9*[-1 1])

%%
FIG = figure (2); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 3 3];
FIG.Name = filename;
movegui(FIG,'center')

HEAD_ALL = cat(2,HEAD{[2:4,6]});
WING_ALL = cat(2,WING{[2:4,6]});

ax = gca; hold on
xlabel('Head (°)')
ylabel('Wing (V)')

xData = HEAD_ALL;
yData = WING_ALL;

[h] = scatplot(xData(:),yData(:));
% delete(h.cb)

% Calculate linear best fit
[r,m,b] = regression(xData',yData','one');
text(ax.XLim(1) + 12, ax.YLim(1) + 0.3,['r =' num2str(r)]);
xFit = linspace(-20,20,4000);
yFit = m*xFit + b;
plot(xFit,yFit,'r','LineWidth',2)

set(ax, 'FontSize', 8, 'LineWidth', 1.5, 'XLim', 20*[-1 1], 'YLim', 9*[-1 1])

%%
FIG = figure (3); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 3 3];
FIG.Name = filename;
movegui(FIG,'center')
for ww = 1:nAmp
    pp = 1;
    for jj = 1
        ax = subplot(1,1,pp) ; hold on        
        ax.XLabel.String = ['Head (' char(176) ')'];
        ax.YLabel.String = 'Wing (V)';
        
       	xData = HEAD{jj,ww};
        yData = WING{jj,ww};
        
        cmap = jet(size(xData,2));

        for kk = 1:size(xData,2)
            h.scat = scatter(xData(:,kk),yData(:,kk),2,'MarkerFaceColor',cmap(kk,:),'MarkerFaceAlpha',0.2,...
                                'MarkerEdgeColor','none');
        end
                
        % Calculate linear best fit
        [r,m,b] = regression(xData,yData,'one');
        text(mean(ax.XLim),mean(ax.YLim),['r =' num2str(r)]);
        xFit = linspace(-20,20,4000);
        yFit = m*xFit + b;
        plot(xFit,yFit,'r','LineWidth',2)
        
        pp =  pp + 1;
    end
end

set(ax, 'FontSize', 8, 'LineWidth', 1.5, 'XLim', 20*[-1 1], 'YLim', 9*[-1 1])

end