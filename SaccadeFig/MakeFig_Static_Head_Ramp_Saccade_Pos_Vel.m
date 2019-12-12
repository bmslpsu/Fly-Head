function [FIG] = MakeFig_Static_Head_Ramp_Saccade_Pos_Vel()
%% MakeFig_Static_Head_Ramp_Saccade_Pos_Vel:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     :   figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

[FILE.Plus,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select positive saccade data', root, 'MultiSelect','off');

[FILE.Minus,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select negative saccade data', root, 'MultiSelect','off');

[FILE.Anti,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select anti saccade data', root, 'MultiSelect','off');

Plus    = load(fullfile(root,FILE.Plus) ,'SACCADE','U','N');
Minus  	= load(fullfile(root,FILE.Minus),'SACCADE','U','N');
Anti    = load(fullfile(root,FILE.Anti) ,'SACCADE','U','N');

CC = repmat(hsv(Anti.N{1,3}/2),2,1);
% Wave = Anti.U{1,3}{1};

%% Combine all anti speeds
FIG = figure (1) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 2*(4/3) 3];
FIG.Name = 'Saccade Velocity';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
clear h

comb_med = cell(Anti.N{1,3}/2,2);
ax = subplot(1,1,1); hold on
for jj = 1:Anti.N{1,3}   
    cent = max(Anti.SACCADE.cIdx{jj},[],'all');
    span = (cent-5):(cent+5);
    
% 	h.trial = plot(1000*Anti.SACCADE.Head.Time{jj},Anti.SACCADE.Head.Velocity{jj}, 'Color', [0.5 0.5 0.5, 0.2]);
%     
% 	PlotPatch(Anti.SACCADE.HeadStats.Velocity(jj).Median(span), Anti.SACCADE.HeadStats.Velocity(jj).STD(span), ...
%         1000*Anti.SACCADE.HeadStats.Time(jj).Median(span), 1, 1, CC(jj,:), [0.7 0.7 0.7], 0.4, 3);
    
    comb_med{jj} = Anti.SACCADE.HeadStats.Velocity(jj).Median(span);

end
set(ax,'FontSize',8,'Color','w','YColor','k','XColor','k','XLim',1000*0.05*[-1 1],'YLim',1100*[-1 1])
xlabel('Time (ms)')
ylabel(['Head Velocity (' char(176) '/s)'])

all_med = cell(1,2);
for jj = 1:2
    all_med{jj} = cat(2,comb_med{:,jj});
end

grand_med = cellfun(@(x) median(x,2), all_med, 'UniformOutput', false);
grand_std = cellfun(@(x) std(x,[],2)   , all_med, 'UniformOutput', false);

for jj = 1:2
    [~,h.anti] = PlotPatch(grand_med{jj}, grand_std{jj}, 1000*Anti.SACCADE.HeadStats.Time(1).Median(span), 1, 1, 'k', [0.7 0.7 0.7], 0.4, 3);
end

for jj = 3 % 30 deg wavelength
    cent = max(Plus.SACCADE.cIdx{jj},[],'all');
    span = (cent-5):(cent+5);
    
    % h.trial = plot(1000*Plus.SACCADE.Head.Time{jj},Plus.SACCADE.Head.Velocity{jj}, 'Color', [0.7*CC(jj,:) , 0.2]);
    
	[~,h.plus] = PlotPatch(Plus.SACCADE.HeadStats.Velocity(jj).Median(span), Plus.SACCADE.HeadStats.Velocity(jj).STD(span), ...
        1000*Plus.SACCADE.HeadStats.Time(jj).Median(span), 1, 1, [0.5 0.5 0.5], [0.7 0.7 0.7], 0.4, 3);
    
 	cent = max(Minus.SACCADE.cIdx{jj},[],'all');
    span = (cent-5):(cent+5);
    
	% h.trial = plot(1000*Minus.SACCADE.Head.Time{jj}, Minus.SACCADE.Head.Velocity{jj}, 'Color', [0.5 0.5 0.5, 0.2]);
    
	[~,h.minus] = PlotPatch(Minus.SACCADE.HeadStats.Velocity(jj).Median(span), Minus.SACCADE.HeadStats.Velocity(jj).STD(span), ...
        1000*Minus.SACCADE.HeadStats.Time(jj).Median(span), 1, 1, [0.5 0.5 0.5], [0.7 0.7 0.7], 0.4, 3);
end

leg = legend([h.anti,h.plus],'Moving','Static');
leg.Box = 'off';

end