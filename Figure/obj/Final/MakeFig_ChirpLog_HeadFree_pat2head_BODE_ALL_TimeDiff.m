function [FIG] = MakeFig_ChirpLog_HeadFree_pat2head_BODE_ALL_TimeDiff()
%% MakeFig_ChirpLog_HeadFree_pat2head_BODE_ALL_TimeDiff:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     : figure handle
%

root = 'H:\DATA\Rigid_Data\';

[CHIRP,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select chirp file', root, 'MultiSelect','on');
CHIRP = cellstr(CHIRP)';

HeadFree = load(fullfile(root,CHIRP{1}),'GRAND','U','N');
clearvars -except HeadFree

%% Head
clearvars -except HeadFree

filename = 'ChirpLog_HeadFree_pat2head_BODE_ALL_new';
catIdx = 7;
xIdx = 1;

freq_bins = ((0:0.4:100)-0.0001)';
GAIN  = nan(length(freq_bins)-1, HeadFree.N.Amp);
PHASE = nan(length(freq_bins)-1, HeadFree.N.Amp);
GSTD  = nan(length(freq_bins)-1, HeadFree.N.Amp);
PSTD  = nan(length(freq_bins)-1, HeadFree.N.Amp);
TD    = nan(length(freq_bins)-1, HeadFree.N.Amp);
TDSTD = nan(length(freq_bins)-1, HeadFree.N.Amp);
for jj = 1:HeadFree.N.Amp % amplitudes
    Fv = HeadFree.GRAND{jj,catIdx}.Mean{2}{1}(:,xIdx); % frequency vector [Hz]
    FREQ_mean = HeadFree.GRAND{jj,catIdx}.Mean{1}{17}(:,xIdx); % complex frequency domain
    FREQ_std = HeadFree.GRAND{jj,catIdx}.Mean{1}{17}(:,xIdx); % complex frequency domain
    
    [newFv,newFREQ_mean,~,~] = DFT_bin(Fv, FREQ_mean, freq_bins);
    [~,newFREQ_std,~,~]      = DFT_bin(Fv, FREQ_std, freq_bins);
    
	GAIN(:,jj)	= abs(newFREQ_mean);
    PHASE(:,jj)	= rad2deg(angle(newFREQ_mean));
    GSTD(:,jj) 	= abs(newFREQ_std);
    PSTD(:,jj) 	= rad2deg(angle(newFREQ_std));
    
    phs_cnd = (PHASE(:,jj) > 0);
    disp(sum(phs_cnd))
   	PHASE(phs_cnd,jj) = PHASE(phs_cnd,jj)  - 360;
    
    TD(:,jj)  	= 1000*( PHASE(:,jj) / 360 ) .* (1./newFv);
    TDSTD(:,jj)	= 1000*( PSTD(:,jj) / 360 ) .* (1./newFv);
end
% indexFv = repmat(1:HeadFree.N.Amp,length(newFv),1);
% phs_cnd = (PHASE > 0) & (repmat(newFv,1,HeadFree.N.Amp) > 5) & ...
%     (( indexFv == 1) | (indexFv == 3) | (indexFv == 2));
% PHASE(phs_cnd) = PHASE(phs_cnd) - 360;

% All Amplitudes
FIG = figure (1); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 HeadFree.N.Amp*2.5 3*2.5];
FIG.Name = filename;
movegui(FIG,'center')
CC = hsv(HeadFree.N.Amp+1);
pp = 1;
ax = gobjects(3,HeadFree.N.Amp);
for jj = 1:HeadFree.N.Amp % amplitudes    
    ax(1,jj) = subplot(3,HeadFree.N.Amp,pp);
        hold on
        ylabel('Gain (�/�)')
        ax(1,jj).YTick = unique(sort([ax(1,jj).YTick ax(1,jj).YLim(2)]));

        h.patch = PlotPatch(GAIN(:,jj), GSTD(:,jj), newFv ,...
                            1, HeadFree.N.fly, CC(jj+1,:), 0.5*CC(jj+1,:), 0.2, 1);
                               
        ax(1,jj).XTick = unique(sort([min(ax(1,jj).XLim) 2:2:12]));
        vel = round(HeadFree.U.Amp{1}(jj)*2*pi*ax(1,jj).XTick);
        velLabel = cellfun(@(x) num2str(x), num2cell(vel), 'UniformOutput', false);
               
    ax(2,jj) = subplot(3,HeadFree.N.Amp,pp + HeadFree.N.Amp);
        hold on
        title([num2str(HeadFree.U.Amp{1}(jj)) , '�'])
        ylabel('Phase Difference (�)')
        % xlabel('Frequency (Hz)')
        
        h.patch = PlotPatch(PHASE(:,jj), PSTD(:,jj), newFv,...
            1, HeadFree.N.fly, CC(jj+1,:), 0.5*CC(jj+1,:), 0.2, 1);
                
        plot([0 12],[0 0],'--k','LineWidth',1)
                
	ax(3,jj) = axes;
        ax(3,jj).Position = ax(1,jj).Position + [0 -0.00 0 0];
        ax(3,jj).Color = 'none';
        ax(3,jj).YAxisLocation = 'right';
        ax(3,jj).YAxis.Color = 'none';
        ax(3,jj).XAxisLocation = 'top';
        ax(3,jj).XTick = ax(1,jj).XTick;
        ax(3,jj).XTickLabels = velLabel;
      	xlabel('Peak Velocity (�/s)')
        
    ax(4,jj) = subplot(3,HeadFree.N.Amp,pp + 2*HeadFree.N.Amp);
        hold on
        ylabel('Time Difference (s)')
        xlabel('Frequency (Hz)')
        ylim([-180 0])
        
        h.patch = PlotPatch(TD(:,jj), TDSTD(:,jj), newFv,...
            1, HeadFree.N.fly, CC(jj+1,:), 0.5*CC(jj+1,:), 0.2, 1);
                
        plot([0 12],[0 0],'--k','LineWidth',1)
        
  	pp = pp + 1;
end
set(ax,'FontSize', 8, 'LineWidth', 1.5, 'XLim',[0.1 12])
set(ax, 'XTick', sort([ax(1,jj).XLim(1),2:2:12]))
set(ax(1,:),'YLim',[0 1])
set(ax(2:3,:), 'YLim', [-300 180], 'YTick', -300:60:180)
linkaxes(ax,'x')

%% All Amplitudes on one
FIG = figure (2); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 1*2.5 5];
FIG.Name = filename;
movegui(FIG,'center')
CC = hsv(HeadFree.N.Amp+1);
pp = 1;
ax = gobjects(2,1);
for jj = 1:HeadFree.N.Amp % amplitudes    
    ax(1) = subplot(2,1,1);
        hold on
        ylabel('Gain (�/�)')
        ax(1,1).YTick = unique(sort([ax(1).YTick ax(1).YLim(2)]));

        h.patch = PlotPatch(GAIN(:,jj), GSTD(:,jj), newFv ,...
                            1, HeadFree.N.fly, CC(jj+1,:), 0.5*CC(jj+1,:), 0.2, 1);
                                              
    ax(2) = subplot(2,1,2);
        hold on
        ylabel('Phase Difference (�)')
        xlabel('Frequency (Hz)')
        
        h.patch = PlotPatch(PHASE(:,jj), PSTD(:,jj), newFv,...
            1, HeadFree.N.fly ,CC(jj+1,:), 0.5*CC(jj+1,:), 0.2, 1);
                
        plot([0 12],[0 0],'--k','LineWidth',1)
        
  	pp = pp + 1;
end
set(ax,'FontSize', 8, 'LineWidth', 1.5, 'XLim',[0.1 12])
set(ax, 'XTick', sort([ax(1).XLim(1),2:2:12]))
set(ax(1,:),'YLim',[0 1])
set(ax(2,:), 'YLim', rad2deg(pi*[-1 1]), 'YTick', -180:60:180)
linkaxes(ax,'x')

end