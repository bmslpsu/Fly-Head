function [FIG] = MakeFig_ChirpLog_HeadFree_pat2wing_BODE_ALL_new3()
%% MakeFig_ChirpLog_HeadFree_pat2wing_BODE_ALL_new3:
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

%% patwing
filename = 'ChirpLog_HeadFree_pat2wing_BODE_ALL_new3';
catIdx = 3;
xIdx = 1;

freq_bins = ((0:0.4:100)-0.0001)';
GAIN  = nan(length(freq_bins)-1, HeadFree.N.Amp);
PHASE = nan(length(freq_bins)-1, HeadFree.N.Amp);
GSTD  = nan(length(freq_bins)-1, HeadFree.N.Amp);
PSTD  = nan(length(freq_bins)-1, HeadFree.N.Amp);
for jj = 1:HeadFree.N.Amp % amplitudes
    Fv = HeadFree.GRAND{jj,catIdx}.Mean{2}{1}(:,xIdx); % frequency vector [Hz]
    
    FREQ_mean = HeadFree.GRAND{jj,catIdx}.Mean{1}{17}(:,xIdx); % complex frequency domain
    FREQ_std = HeadFree.GRAND{jj,catIdx}.Mean{1}{17}(:,xIdx); % complex frequency domain
    
    gain_mean  = HeadFree.GRAND{jj,catIdx}.Mean{1}{2}(:,xIdx); % gain
    gain_std  = HeadFree.GRAND{jj,catIdx}.Mean{1}{2}(:,xIdx); % gain
    phase_mean = HeadFree.GRAND{jj,catIdx}.CircMean{7}{2}(:,xIdx); % phase
    phase_std = HeadFree.GRAND{jj,catIdx}.CircSTD{7}{2}(:,xIdx); % phase
    
    [newFv,newFREQ_mean,~,~] = DFT_bin(Fv, FREQ_mean, freq_bins);
    [~,newFREQ_std,~,~]      = DFT_bin(Fv, FREQ_std, freq_bins);
    
 	[~,newGain_mean,~,~] = DFT_bin(Fv, gain_mean, freq_bins);
    [~,newGain_std,~,~]  = DFT_bin(Fv, gain_std, freq_bins);
	[~,newPhase_mean,~,~] = DFT_bin(Fv, gain_mean, freq_bins);
    [~,newPhase_std,~,~]  = DFT_bin(Fv, phase_std, freq_bins);
    
	GAIN(:,jj)	= newGain_mean;
    PHASE(:,jj)	= rad2deg(angle(newFREQ_mean));
    GSTD(:,jj)	= newGain_std;
    PSTD(:,jj)	= rad2deg(angle(newFREQ_std));
    
    GAIN(1,jj)	= GAIN(2,jj);

    % GAIN(:,jj)     = abs(newFREQ_mean);
    % PHASE(:,jj)    = rad2deg(angle(newFREQ_mean));
    % GSTD(:,jj)     = abs(newFREQ_std);
    % PSTD(:,jj)     = rad2deg(angle(newFREQ_std));
end
indexFv = repmat(1:HeadFree.N.Amp,length(newFv),1);
phs_cnd = (PHASE > 20) & (repmat(newFv,1,HeadFree.N.Amp) > 2) & ...
    (( indexFv == 1) | (indexFv == 4) | (indexFv == 2) | (indexFv == 3));
PHASE(phs_cnd) = PHASE(phs_cnd) - 360;

phs_cnd_2 = (PHASE > 20) & (repmat(newFv,1,HeadFree.N.Amp) > 1.5) & (indexFv == 2);
PHASE(phs_cnd_2) = PHASE(phs_cnd_2) - 360;

%% All Amplitudes
FIG = figure (1); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 HeadFree.N.Amp*2.5 5];
FIG.Name = filename;
movegui(FIG,'center')
CC = hsv(HeadFree.N.Amp+1);
pp = 1;
ax = gobjects(3,HeadFree.N.Amp);
for jj = 1:HeadFree.N.Amp % amplitudes    
    ax(1,jj) = subplot(2,HeadFree.N.Amp,pp);
        hold on
        ylabel('Gain (°/°)')
        ax(1,jj).YTick = unique(sort([ax(1,jj).YTick ax(1,jj).YLim(2)]));

        h.patch = PlotPatch(GAIN(:,jj), GSTD(:,jj), newFv ,...
                    1, HeadFree.N.fly, CC(jj+1,:), 0.5*CC(jj+1,:), 0.2, 1);
                               
        ax(1,jj).XTick = unique(sort([min(ax(1,jj).XLim) 2:2:12]));
        vel = round(HeadFree.U.Amp{1}(jj)*2*pi*ax(1,jj).XTick);
        velLabel = cellfun(@(x) num2str(x), num2cell(vel), 'UniformOutput', false);
               
    ax(2,jj) = subplot(2,HeadFree.N.Amp,pp + HeadFree.N.Amp);
        hold on
        title([num2str(HeadFree.U.Amp{1}(jj)) , '°'])
        ylabel('Phase Difference (°)')
        xlabel('Frequency (Hz)')
        
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
      	xlabel('Peak Velocity (°/s)')
        
  	pp = pp + 1;
end
set(ax,'FontSize', 8, 'LineWidth', 1.5, 'XLim',[0.1 12])
set(ax, 'XTick', sort([ax(1,jj).XLim(1),2:2:12]))
set(ax(1,:),'YLim',[0 0.5])
set(ax(1,:),'YTick',0:0.1:0.5)
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
        ylabel('Gain (°/°)')
        ax(1,1).YTick = unique(sort([ax(1).YTick ax(1).YLim(2)]));

        h.patch = PlotPatch(GAIN(:,jj), GSTD(:,jj), newFv ,...
                            1, HeadFree.N.fly, CC(jj+1,:), 0.5*CC(jj+1,:), 0.2, 1);
                                              
    ax(2) = subplot(2,1,2);
        hold on
        ylabel('Phase Difference (°)')
        xlabel('Frequency (Hz)')
        
        h.patch = PlotPatch(PHASE(:,jj), PSTD(:,jj), newFv,...
            1, HeadFree.N.fly ,CC(jj+1,:), 0.5*CC(jj+1,:), 0.2,1);
                
        plot([0 12],[0 0],'--k','LineWidth',1)
        
  	pp = pp + 1;
end
set(ax,'FontSize', 8, 'LineWidth', 1.5, 'XLim',[0.1 12])
set(ax, 'XTick', sort([ax(1,1).XLim(1),2:2:12]))
set(ax(1,:),'YLim',[0 0.5])
set(ax(1,:),'YTick',0:0.1:0.5)
set(ax(2,:), 'YLim', [-300 180], 'YTick', -300:60:180)
linkaxes(ax,'x')

end