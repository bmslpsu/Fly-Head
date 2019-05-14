function [FIG] = MakeFig_Ramp_HeadFree_Wing_Freq()
%% MakeFig_Ramp_HeadFree_Wing_Freq: wing frequency spectrum plot for ramp
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'F:\DATA\Rigid_Data';

spatFreq = [22.5 30 60 0];
nFreq = length(spatFreq);
HeadFree = cell(nFreq,1);
for ww = 1:nFreq
    HeadFree{ww} = load(fullfile(root, ['Ramp_HeadFree_' num2str(spatFreq(ww)) '_DATA.mat']),...
        'TRIAL','FLY','GRAND','U','N'); % load data structure
end

filename = 'Ramp_HeadFree_Wing_Freq'; % name of figure to save
figNum = 1;
catIdx = 3;
xIdx = 1;

FIG = figure (figNum); % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1100 800];
FIG.Name = filename;
hold on

subIdx = reshape(1:(nFreq*5), 4, 5);

% Trials
pp = 1;
for ww = 1:nFreq % spatial frequencies
    for jj = 1:HeadFree{ww}.N{1,3} % speeds
        subplot(nFreq,HeadFree{ww}.N{1,3},subIdx(pp)) ; hold on ; xlim([0 30]) ; ylim([0 1])
        for kk = 1:HeadFree{ww}.N{1,1} % flys
            for ii = 1: size(HeadFree{ww}.TRIAL{kk,jj},1) % trials
                h.Trial = plot(HeadFree{ww}.TRIAL{kk,jj}{ii,catIdx}.Fv(:,xIdx),HeadFree{ww}.TRIAL{kk,jj}{ii,catIdx}.Mag(:,xIdx),...
                    '-','Color',[0.5 0.5 0.5 0.25],'LineWidth',1);
            end
        end
        pp = pp + 1;
    end
end

% Fly Stats
pp = 1;
for ww = 1:nFreq % spatial frequencies
    for jj = 1:HeadFree{ww}.N{1,3} % speeds
        subplot(nFreq,HeadFree{ww}.N{1,3},subIdx(pp)) ; hold on ; xlim([0 30]) ; ylim([0 1])
        for kk = 1:HeadFree{ww}.N{1,1} % flys
            h.Fly = plot(HeadFree{ww}.FLY{jj}{kk,catIdx}.Mean{7}(:,xIdx),HeadFree{ww}.FLY{jj}{kk,catIdx}.Mean{8}(:,xIdx),'-','LineWidth',2);
            h.Fly.Color(4) = 0.5;
        end
        pp = pp + 1;
    end
end


% Grand Stats
pp = 1;
for ww = 1:nFreq % spatial frequencies
    freq = spatFreq(ww);
    for jj = 1:HeadFree{ww}.N{1,3} % speeds
     	subplot(nFreq,HeadFree{ww}.N{1,3},subIdx(pp)) ; hold on ; xlim([0 40]) ; ylim([0 0.2])
%         grid on ; grid minor
        vel = HeadFree{1}.U{1,3}{1}(jj);
        tmpFreq = vel/freq;
        if pp<=5
            title([num2str(vel) char(176) '/s'],'FontSize',15,'FontWeight','bold')
        end
        if mod(pp-1,5)==0
            ylabel([num2str(freq) char(176)],'FontSize',15,'FontWeight','bold')
        else
            yticks('')
        end
        if pp>15
            xlabel('Frequency (Hz)','FontSize',12,'FontWeight','normal')
        else
            xticks('')
        end
        
        text(15,0.1,[num2str(tmpFreq) ' cyc/s'])
        
        h.patch = PlotPatch(HeadFree{ww}.GRAND{jj,catIdx}.Mean{2}{8}(:,xIdx),HeadFree{ww}.GRAND{jj,catIdx}.STD{2}{8}(:,xIdx),...
            HeadFree{ww}.GRAND{jj,catIdx}.Mean{2}{7}(:,xIdx),2,HeadFree{ww}.N{1,1},'k',[0.4 0.4 0.6],0.5,3);
        
        pp = pp + 1;
    end
end

saveas(FIG,['F:\DATA\Rigid_Data\FIGURE\' filename '.fig']); % save .fig file
% print(FIG,['F:\DATA\FIGURE\' filename],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
end