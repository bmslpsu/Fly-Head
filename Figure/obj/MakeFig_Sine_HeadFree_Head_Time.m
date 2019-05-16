function [] = MakeFig_Sine_HeadFree_Head_Time()
%% MakeFig_Sine_HeadFree_Head_Time:
%   INPUTS:
%       -
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

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
    HeadFree{ww} = load(fullfile(root,FILES{ww}),'TRIAL','FLY','GRAND','U','N');
end

filename = 'Sine_HeadFree_Head_Time'; % name of figure to save

hold on
catIdx = 2;
xIdx = 1;
figNum = 1;

FIG = figure (figNum); clf % figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1100 800];
movegui(FIG,'center')
FIG.Name = filename;
for ww = 1:nAmp
   FIG.Name = [FIG.Name '_' num2str(Amp(ww))];  
end

% Trials
for ww = 1:nAmp % amplitudes
    pp = 1;
    for jj = 1:HeadFree{ww}.N{1,3} % frequencies
        subplot(HeadFree{ww}.N{1,3},1,pp) ; hold on ; xlim([0 10]) ; ylim(20*[-1 1])
        for kk = 1:HeadFree{ww}.N{1,1} % flys
            for ii = 1:size(HeadFree{ww}.TRIAL{kk,jj},1) % trials
                h.Trial = plot(HeadFree{ww}.TRIAL{kk,jj}{ii,catIdx}.Time,HeadFree{ww}.TRIAL{kk,jj}{ii,catIdx}.X(:,xIdx),...
                    '-','Color',[0.5 0.5 0.5 0.25],'LineWidth',1);
            end
        end
        pp = pp + 1;
    end
end

% Fly Stats
for ww = 1:nAmp % amplitudes
    pp = 1;
    for jj = 1:HeadFree{ww}.N{1,3} % frequencies
        subplot(HeadFree{ww}.N{1,3},1,pp) ; hold on ; xlim([0 10]) ; ylim(20*[-1 1])
        for kk = 1:HeadFree{ww}.N{1,1} % flys
            h.Fly = plot(HeadFree{ww}.FLY{jj}{kk,catIdx}.Mean{5},HeadFree{ww}.FLY{jj}{kk,catIdx}.Mean{6}(:,xIdx),'-','LineWidth',2);
            h.Fly.Color(4) = 0.5;
        end
        pp = pp + 1;
    end
end

% Grand Stats
for ww = 1:nAmp % amplitudes 
    pp = 1;
    for jj = 1:HeadFree{ww}.N{1,3} % frequencies
        freq = HeadFree{ww}.U{1,3}{1}(jj);
     	subplot(HeadFree{ww}.N{1,3},1,pp)  ; hold on ; xlim([0 10]) ; ylim(20*[-1 1])

     	title([num2str(freq) ' Hz'],'FontSize',15,'FontWeight','bold')
        
       	ylabel(['Head ($^{\circ}$)'],'interpreter','latex','FontSize',15)
        
        if pp==HeadFree{ww}.N{1,3}
            xlabel('Time (s)','interpreter','latex','FontSize',15)
        else
            xticks('')
        end
               
        h.patch = PlotPatch(HeadFree{ww}.GRAND{jj,catIdx}.Mean{2}{6}(:,xIdx),HeadFree{ww}.GRAND{jj,catIdx}.STD{2}{6}(:,xIdx),...
            HeadFree{ww}.GRAND{jj,catIdx}.Mean{2}{5},2,HeadFree{ww}.N{1,1},'k',[0.4 0.4 0.6],0.5,3);
        
        pp = pp + 1;
    end
end

figdir = 'H:\DATA\Rigid_Data\FIGURE\';
saveas(FIG,[figdir FIG.Name '.fig']); % save .fig file
%print (FIG,[figdir FIG.Name],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
end