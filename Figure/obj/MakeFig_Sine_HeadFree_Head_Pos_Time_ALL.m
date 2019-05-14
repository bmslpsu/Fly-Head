function [FIG] = MakeFig_Sine_HeadFree_Head_Pos_Time_ALL()
%% MakeFig_Sine_HeadFree_Head_Pos_Time_ALL:
%   INPUTS:
%       -
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
    HeadFree{ww} = load(fullfile(root,FILES{ww}),'TRIAL','FLY','GRAND','U','N');
end

figNum = 1;
filename = 'Sine_HeadFree_Head_Pos_Time_ALL'; % name of figure to save

FIG = figure (figNum); clf% figure handle
FIG.Color = 'w';
FIG.Position = [100 100 1100 800];
FIG.Name = filename;
for ww = 1:nAmp
   FIG.Name = [FIG.Name '_' num2str(Amp(ww))];  
end
hold on

catIdx = 2;
xIdx = 1;

CLR         = jet(nAmp);
legList     = cellstr(num2str(Amp))';
hh          = cell(nAmp,1);
% Grand Stats
for ww = 1:nAmp % amplitudes 
    pp = 1;
    for jj = 1:HeadFree{ww}.N{1,3} % frequencies
        freq = HeadFree{ww}.U{1,3}{1}(jj);
     	subplot(HeadFree{ww}.N{1,3},1,pp)  ; hold on ; xlim([0 10]) ; ylim(20*[-1 1])

     	title([num2str(freq) ' Hz'],'Interpreter','latex','FontSize',15)
        
       	ylabel('Head ($^{\circ}$)','Interpreter','latex','FontSize',12)
        
        if pp==HeadFree{ww}.N{1,3}
            xlabel('Time (s)','Interpreter','latex','FontSize',15)
        else
            xticks('')
        end
               
       [~,h.patch] = PlotPatch(HeadFree{ww}.GRAND{jj,catIdx}.Mean{2}{6}(:,xIdx),HeadFree{ww}.GRAND{jj,catIdx}.STD{2}{6}(:,xIdx),...
            HeadFree{ww}.GRAND{jj,catIdx}.Mean{2}{5},2,100,CLR(ww,:),[0.4 0.4 0.6],0.5,2);
        hh{ww} = h.patch;
                
        pp = pp + 1;
    end
end
subplot(HeadFree{ww}.N{1,3},1,1)
legend(cat(1,hh{:}),legList)

figdir = 'F:\DATA\Rigid_Data\FIGURE\';
saveas(FIG,[figdir FIG.Name '.fig']); % save .fig file
% print (FIG,[figdir FIG.Name],'-dpdf','-r600','-bestfit') % save as publication quality .pdf
end