function [FIG] = MakeFig_Sine_HeadFree_patvshead_TimeScatter_ALL()
%% MakeFig_Sine_HeadFree_headvswing_TimeScatter_ALL:
%   INPUTS:
%      -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'F:\DATA\Rigid_Data\';

% Select files
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
    HeadFree{ww} = load(fullfile(root,FILES{ww}),'TRIAL','U','N');
end

figNum = 1;
filename = 'Sine_HeadFree_patvshead_TimeScatter_ALL'; % name of figure to save
catIdx = [1 2];
xIdx = 1;

% Store data by amplitude
PAT = cell(nAmp,HeadFree{1}.N{1,3});
HEAD = PAT;
pp = 1;
for ww = 1:nAmp
    for jj = 1:HeadFree{ww}.N{1,3}
        for kk = 1:HeadFree{ww}.N{1,1}
            for ii = 1:size(HeadFree{ww}.TRIAL{kk,jj},1)
                PAT{ww,jj}(:,end+1) = HeadFree{ww}.TRIAL{kk,jj}{ii,catIdx(1)}.X(:,xIdx);
                HEAD{ww,jj}(:,end+1) = HeadFree{ww}.TRIAL{kk,jj}{ii,catIdx(2)}.X(:,xIdx);
                pp = pp + 1; 
            end
        end
    end
end

%
FIG = figure (figNum); clf % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1200 220*nAmp];
FIG.Name = filename;
movegui(FIG,'center')
hold on

pp = 1;
for ww = 1:nAmp
    for jj = 1:HeadFree{ww}.N{1,3}
        ax = subplot(nAmp,HeadFree{ww}.N{1,3},pp) ; hold on
        xlim(20*[-1 1])
        ylim(8*[-1 1])
        if pp<=HeadFree{ww}.N{1,3}
            title([num2str(HeadFree{ww}.U{1,3}{1}(jj)) ' Hz'],'Interpreter','latex','FontSize',15)
        end
            
        if pp>(nAmp-1)*HeadFree{ww}.N{1,3}
            xlabel('Pattern ($^{\circ}$)','Interpreter','latex','FontSize',12)
        else
            xticks('')
        end
        
        if round((pp-1)/HeadFree{ww}.N{1,3})==((pp-1)/HeadFree{ww}.N{1,3})
           ylabel({[num2str(Amp(ww)) '$^{\circ}$'],'Head ($^{\circ}$)'},'Interpreter','latex','FontSize',12)
        else
            yticks('')
        end

        xData = PAT{ww,jj}(:);
        yData = HEAD{ww,jj}(:);

        [h] = scatplot(xData,yData);
        delete(h.cb)
%         scatter(xData(:,1),yData(:,1))
        
        % Calculate linear best fit
        [r,m,b] = regression(xData,yData,'one');
        text(0,-7,['r =' num2str(r)]);
        xFit = linspace(-20,20,4000);
        yFit = m*xFit + b;
        plot(xFit,yFit,'r','LineWidth',5)
        
    pp = pp + 1;
%     lsline
    end
end

figdir = 'F:\DATA\Rigid_Data\FIGURE\';
saveas(FIG,[figdir FIG.Name '.fig']); % save .fig file
% print (FIG,[figdir FIG.Name],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
end