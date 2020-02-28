function [] = MakeFig_SOS_Correlation()
%% MakeFig_SOS_Correlation:

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

stimIdx = 1;
headIdx = 2;
wingIdx = 3;
errIdx  = 4;
xIdx = 1;

nFreq = HeadFree{1}.N{1,3};

STIM = cell(nFreq,nAmp);
HEAD = cell(nFreq,nAmp);
WING = cell(nFreq,nAmp);
ERR  = cell(nFreq,nAmp);
for ww = 1:nAmp
    nFly = HeadFree{ww}.N{1,1};
    for jj = 1:nFreq
        pp = 1;
        for kk = 1:nFly
            for ii = 1:size(HeadFree{ww}.TRIAL{kk,jj},1) % trials
                STIM{jj,ww}(:,pp) = HeadFree{ww}.TRIAL{kk,jj}{ii,stimIdx}.X(:,xIdx);
                HEAD{jj,ww}(:,pp) = HeadFree{ww}.TRIAL{kk,jj}{ii,headIdx}.X(:,xIdx);
                WING{jj,ww}(:,pp) = HeadFree{ww}.TRIAL{kk,jj}{ii,wingIdx}.X(:,xIdx);
                ERR{jj,ww}(:,pp)  = HeadFree{ww}.TRIAL{kk,jj}{ii,errIdx }.X(:,xIdx);
                
%                 cla
%                 time = HeadFree{ww}.TRIAL{kk,jj}{ii,stimIdx}.Time;
%                 head = HEAD{jj,ww}(:,pp);
%                 wing = 5*WING{jj,ww}(:,pp);
%                 stim = STIM{jj,ww}(:,pp);
%                 hold on
%                 plot(time, stim,'k')
%                 plot(time, head,'b')
%                 plot(time, wing,'r')
%                 
%                 pause
                
                pp =  pp + 1;
            end
        end
    end
end

%% Head vs Wing
fig = figure (1); clf
set(fig,'Color','w','Units','inches', 'Position', [2 2 3 3], 'Name', 'Head vs Wing');
movegui(fig,'center')

ax = gobjects(nAmp,nFreq);
pp = 1;
for ww = 1:nAmp
    for jj = 1:nFreq
     	xData = HEAD{jj,ww}(:);
        yData = WING{jj,ww}(:);
        
    	% Calculate linear best fit
        [r,m,b] = regression(xData,yData,'one');
        xFit = linspace(-20,20,4000);
        yFit = m*xFit + b;
        
        ax(pp) = subplot(nAmp,nFreq,pp) ; hold on        
            xlabel({'Head (°)'})
            ylabel('Wing (V)');
%             [h] = scatplot(xData,yData,[],[],[],[],1,[]);
%             delete(h.cb)

            H = densityplot(xData,yData);
            
        % plot(xFit,yFit,'r','LineWidth',2)
        % title(['r = ' num2str(r)])
        
        pp =  pp + 1;
    end
end
% set(ax, 'FontSize', 8, 'LineWidth', 1.5, 'XLim', 20*[-1 1], 'YLim', 9*[-1 1])


%% Stim vs Head
fig = figure (2); clf
set(fig,'Color','w','Units','inches', 'Position', [2 2 3 3], 'Name', 'Stimulus vs Head');
movegui(fig,'center')

ax = gobjects(nAmp,nFreq);
pp = 1;
for ww = 1:nAmp
    for jj = 1:nFreq
     	xData = STIM{jj,ww}(:);
        yData = HEAD{jj,ww}(:);
        
    	% Calculate linear best fit
        [r,m,b] = regression(xData,yData,'one');
        xFit = linspace(-20,20,4000);
        yFit = m*xFit + b;
        
        ax(pp) = subplot(nAmp,nFreq,pp) ; hold on        
            xlabel('Stimulus (°)')
            ylabel('Head (°)');
%             [h] = scatplot(xData,yData);
%             delete(h.cb)

        H = densityplot(xData,yData);

        plot(xFit,yFit,'r','LineWidth',2)
        title(['r = ' num2str(r)])
        
        pp =  pp + 1;
    end
end
% set(ax, 'FontSize', 8, 'LineWidth', 1.5, 'XLim', 20*[-1 1], 'YLim', 20*[-1 1])

%% Stim vs Wing
fig = figure (3); clf
set(fig,'Color','w','Units','inches', 'Position', [2 2 3 3], 'Name', 'Stimulus vs Wing');
movegui(fig,'center')

ax = gobjects(nAmp,nFreq);
pp = 1;
for ww = 1:nAmp
    for jj = 1:nFreq
     	xData = STIM{jj,ww}(:);
        yData = WING{jj,ww}(:);
        
    	% Calculate linear best fit
        [r,m,b] = regression(xData,yData,'one');
        xFit = linspace(-20,20,4000);
        yFit = m*xFit + b;
        
        ax(pp) = subplot(nAmp,nFreq,pp) ; hold on        
            xlabel('Stimulus (°)')
            ylabel('Wing (V)');
%             [h] = scatplot(xData,yData);
%             delete(h.cb)
            
            H = densityplot(xData,yData);

        plot(xFit,yFit,'r','LineWidth',2)
        title(['r = ' num2str(r)])
        
        pp =  pp + 1;
    end
end
set(ax, 'FontSize', 8, 'LineWidth', 1.5, 'XLim', 20*[-1 1], 'YLim', 9*[-1 1])

end