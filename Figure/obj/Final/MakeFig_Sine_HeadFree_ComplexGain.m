function [FIG] = MakeFig_Sine_HeadFree_ComplexGain()
%% MakeFig_Sine_HeadFree_ComplexGain:
%   INPUTS:
%       -
%   OUTPUTS:
%       -
%
root = 'H:\DATA\Rigid_Data\';

[FILES,~] = uigetfile({'*.mat', 'DAQ-files'},...
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
    HeadFree{ww} = load(fullfile(root,FILES{ww}),'TRIAL','GRAND','U','N','T');
end

%% Complex Gain Calculations
clearvars -except nAmp Amp HeadFree

filename = 'Sine_HeadFree_ComplexGain';

catIdx = 5;
xIdx = 1;

Freq = HeadFree{1}.U{1,3}{1}';
nFreq = HeadFree{1}.N{1,3};
fLabel = cellfun(@(x) [num2str(x) ' Hz  : '], num2cell(Freq),'UniformOutput',false);

vel = nan(1,nFreq);
velIdx = [4,4,3,2,1,1];
Vel = repmat(2*pi*HeadFree{1}.U{1,3}{1}',nAmp,1);
for ww = 1:nAmp
   Vel(ww,:) = round(Amp(ww)*Vel(ww,:));
end
for jj = 1:nFreq
    vel(jj) = Vel(velIdx(jj),jj);
end
vLabel = cellfun(@(x) strcat(num2str(x),[' ' char(176) '/s']), num2cell(Vel),'UniformOutput',false);
vLabel_Norm = cellfun(@(x) strcat(num2str(x),[' ' char(176) '/s']), num2cell(vel),'UniformOutput',false);

legLabel = cellfun(@(x,y) [x ' ' y], repmat(fLabel,nAmp,1), vLabel,'UniformOutput',false);
legLabel_Norm = cellfun(@(x,y) [x ' ' y], fLabel, vLabel_Norm,'UniformOutput',false);

cList = prism(nFreq);

% Get real & imaginary parts
maxTrial = nan(nAmp,nFreq);
Real = cell(nAmp,1);
Imag = Real;
for ww = 1:nAmp
    Real{ww} = nan(100,nFreq);
  	Imag{ww} = nan(100,nFreq);
    for jj = 1:nFreq
        pp = 1;
        for kk = 1:HeadFree{ww}.N{1,1}
            for ii = 1:size(HeadFree{ww}.TRIAL{kk,jj},1)
                R = real(HeadFree{ww}.TRIAL{kk,jj}{ii,catIdx}.IOCmplxGain(:,xIdx));
                I = imag(HeadFree{ww}.TRIAL{kk,jj}{ii,catIdx}.IOCmplxGain(:,xIdx));
                
                if jj==1 && I<-0.05
                    I = -I;
                    R = -R;
                elseif jj==5 && (R<0)
                 	I = -I;
                    R = -R;
                end
                
                Real{ww}(pp,jj) = R;
                Imag{ww}(pp,jj) = I;
                
                pp = pp + 1;
            end
        end
        maxTrial(ww,jj) = pp - 1;
    end
    maxFreq = max(maxTrial(ww,:));
    Real{ww} = Real{ww}(1:maxFreq,:);
    Imag{ww} = Imag{ww}(1:maxFreq,:);
end
maxALL = max(max(maxTrial));

CmplxGain = cellfun(@(x,y) (x + 1i*y), Real, Imag, 'UniformOutput', false);
Gain = cellfun(@(x) abs(x), CmplxGain, 'UniformOutput', false);
Phase = cellfun(@(x) rad2deg(angle(x)), CmplxGain, 'UniformOutput', false);

for ww = 1:nAmp
   for jj = 1:nFreq
       phase = Phase{ww}(:,jj);
       if (jj==5) && (Amp(ww)==18.75)
           cnd = phase>0;
           Phase{ww}(cnd,jj) = phase(cnd) - 180;
       end
   end    
end

GAIN = cellfun(@(x) nanmean(x,1), Gain, 'UniformOutput', false);
GAIN_STD = cellfun(@(x) nanstd(x,[],1), Gain, 'UniformOutput', false);
PHASE = cell(nAmp,1);
PHASE_STD = cell(nAmp,1);
for ww = 1:nAmp
   for jj = 1:nFreq
       phase = deg2rad(Phase{ww}(:,jj));
       [PHASE{ww}(:,jj)] = rad2deg(circ_mean(phase(~isnan(phase)),[],1));
       [~,temp] = circ_std(phase(~isnan(phase)),[],[],1);
       PHASE_STD{ww}(:,jj) = rad2deg(temp);
   end    
end

GAIN = cat(1,GAIN{:});
PHASE = cat(1,PHASE{:});
GAIN_STD = cat(1,GAIN_STD{:});
PHASE_STD = cat(1,PHASE_STD{:});

Real_Norm = nan(maxALL,nFreq);
Imag_Norm = Real_Norm;
for ww = 1:nAmp
    for jj = 1:nFreq
        if ww==velIdx(jj)
            Real_Norm(1:maxTrial(ww,jj),jj) = Real{ww}(1:maxTrial(ww,jj),jj);
            Imag_Norm(1:maxTrial(ww,jj),jj) = Imag{ww}(1:maxTrial(ww,jj),jj);
        end
    end
end

REAL = cellfun(@(x) nanmean(x,1), Real, 'UniformOutput', false);
IMAG = cellfun(@(x) nanmean(x,1), Imag, 'UniformOutput', false);
REAL = cat(1,REAL{:});
IMAG = cat(1,IMAG{:});

% REAL_STD = cellfun(@(x) nanstd(x,[],1), Real, 'UniformOutput', false);
% IMAG_STD = cellfun(@(x) nanstd(x,[],1), Imag, 'UniformOutput', false);
% REAL_STD = cat(1,REAL_STD{:});
% IMAG_STD = cat(1,IMAG_STD{:});

REAL_NORM = nan(1,nFreq);
IMAG_NORM = REAL_NORM;
for jj = 1:nFreq
    REAL_NORM(jj) = REAL(velIdx(jj),jj);
    IMAG_NORM(jj) = IMAG(velIdx(jj),jj);
end

CmplxGain_Norm      = Real_Norm + 1i*Imag_Norm;
Gain_Norm           = abs(CmplxGain_Norm);
Phase_Norm          = rad2deg(angle(CmplxGain_Norm));
GAIN_NORM           = nanmean(Gain_Norm,1);
PHASE_NORM          = nanmean(Phase_Norm,1);
GAIN_NORM_STD       = nanstd(Gain_Norm,[],1);
PHASE_NORM_STD      = nanstd(Phase_Norm,[],1);

% gains = 0.05:0.05:0.1;
gains = 0.2:0.2:1;

%% Complex Gain: one amplitude
FIG = figure (10); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 4 4];
movegui(FIG,'center')
FIG.Name = filename;
amp = 1;

[ax,h] = ComplexAxes(gains);
% h.circle(end).Color = 'r';
% h.origin(1).Color = 'r';
ax.Title.String = [num2str(Amp(amp)) char(176)];

for jj = 1:nFreq
    h.trial = scatter(Real{amp}(:,jj), Imag{amp}(:,jj), 40, 'o','MarkerEdgeColor','k',...
        'MarkerFaceColor',cList(jj,:), 'MarkerFaceAlpha', 0.65, 'LineWidth', 0.5);
end

for jj = 1:nFreq
	h.r = plot([0 REAL(amp,jj)],[0 IMAG(amp,jj)],'Color',[0 0 0 1],'LineWidth',1);

    rSTD = PolarSTD(Real{amp}(:,jj),Imag{amp}(:,jj),[REAL(amp,jj) IMAG(amp,jj)]);

    [h.std] = draw_ellipse([REAL(amp,jj) IMAG(amp,jj)], 2*rSTD, 0.5, 0, 90, cList(jj,:)); hold on
    h.std{1}.FaceAlpha = 0.5;
    h.std{2}.LineWidth = 1;
    for kk = 3:length(h.std)
       delete(h.std{kk}) 
    end

    h.grand = scatter(REAL(amp,jj),IMAG(amp,jj),1,'o','MarkerEdgeColor','k','MarkerFaceColor',cList(jj,:),...
        'MarkerFaceAlpha',1,'LineWidth',1.5);

    h.leg = scatter(REAL(amp,jj),IMAG(amp,jj),10,'o','MarkerEdgeColor','k','MarkerFaceColor','k',...
        'MarkerFaceAlpha',1,'LineWidth',1.5);

    hh(jj) = h.grand;

    uistack(h.grand,'top')
end

%     scatter(1,0,40,'o','MarkerEdgeColor','r','MarkerFaceColor','r',...
%         'MarkerFaceAlpha',1,'LineWidth',1.5);

leg = legend(hh,legLabel{amp,:});
leg.Title.String = 'Frequency / Velocity';
leg.Location = 'northwest';
leg.Position = leg.Position + [-0.2 0 0 0];
legend boxoff

%% Complex Gain: Normalized Amplitude
FIG = figure (1); clf
FIG.Color = 'w';
FIG.Position = [100 100 700 700];
movegui(FIG,'center')
FIG.Name = filename;
for ww = 1:nAmp
   FIG.Name = [FIG.Name '_' num2str(Amp(ww))];  
end

[~,h.ax] = ComplexAxes(gains);
h.ax.origin(1).Color = 'r';
h.ax.circle(end).Color = 'r';

for jj = 1:nFreq
    h.trial = scatter(Real_Norm(:,jj), Imag_Norm(:,jj), 40, 'o','MarkerEdgeColor','k',...
        'MarkerFaceColor',cList(jj,:), 'MarkerFaceAlpha', 0.65, 'LineWidth', 0.5);
    
% 	h.rr = plot([0 REAL_NORM(jj)],[0 IMAG_NORM(jj)],'Color',[0 0 0 1],'LineWidth',1);
    
    rSTD = PolarSTD(Real_Norm(:,jj),Imag_Norm(:,jj),[REAL_NORM(jj) IMAG_NORM(jj)]);
    
	[h.std] = draw_ellipse([REAL_NORM(jj) IMAG_NORM(jj)], 3*rSTD, 0.5, 0, 90, cList(jj,:)); hold on
    h.std{1}.FaceAlpha = 0.2;
    for kk = 2:length(h.std)
       delete(h.std{kk}) 
    end
    
    h.error = plot([1 REAL_NORM(jj)],[0 IMAG_NORM(jj)],'Color','k','LineWidth',1);
    
    h.grand = scatter(REAL_NORM(jj),IMAG_NORM(jj),1,'o','MarkerEdgeColor','k','MarkerFaceColor',cList(jj,:),...
        'MarkerFaceAlpha',1,'LineWidth',1.5);
    
    h.leg = scatter(REAL_NORM(jj),IMAG_NORM(jj),40,'o','MarkerEdgeColor','k','MarkerFaceColor','k',...
        'MarkerFaceAlpha',1,'LineWidth',1.5);

    hh(jj) = h.grand;
end

scatter(1,0,100,'o','MarkerEdgeColor','r','MarkerFaceColor','r',...
        'MarkerFaceAlpha',1,'LineWidth',1.5);
text(1.04,0.05,'1 + 0i')

leg = legend(hh,legLabel_Norm);
leg.Title.String = 'Frequency / Velocity';
leg.Location = 'northwest';
legend boxoff

%% Complex Gain: All Amplitude
FIG = figure (2); clf
FIG.Color = 'w';
FIG.Position = [100 100 1400 800];
movegui(FIG,'center')
FIG.Name = filename;
for ww = 1:nAmp
   FIG.Name = [FIG.Name '_' num2str(Amp(ww))];  
end

for ww = 1:nAmp
    ax = subplot(ceil(nAmp/2),2,ww);
    [ax,h] = ComplexAxes(gains);
    h.circle(end).Color = 'r';
    h.origin(1).Color = 'r';
    ax.Title.String = [num2str(Amp(ww)) char(176)];

  	for jj = 1:nFreq
        h.trial = scatter(Real{ww}(:,jj), Imag{ww}(:,jj), 40, 'o','MarkerEdgeColor','k',...
            'MarkerFaceColor',cList(jj,:), 'MarkerFaceAlpha', 0.65, 'LineWidth', 0.5);
    end
    
    for jj = 1:nFreq
%         h.r = plot([0 REAL(ww,jj)],[0 IMAG(ww,jj)],'Color',[0 0 0 1],'LineWidth',1);
        
        rSTD = PolarSTD(Real{ww}(:,jj),Imag{ww}(:,jj),[REAL(ww,jj) IMAG(ww,jj)]);
        
     	[h.std] = draw_ellipse([REAL(ww,jj) IMAG(ww,jj)], 3*rSTD, 0.5, 0, 90, cList(jj,:)); hold on
        h.std{1}.FaceAlpha = 0.2;
        for kk = 2:length(h.std)
           delete(h.std{kk}) 
        end
        
        h.grand = scatter(REAL(ww,jj),IMAG(ww,jj),1,'o','MarkerEdgeColor','k','MarkerFaceColor',cList(jj,:),...
            'MarkerFaceAlpha',1,'LineWidth',1.5);

        h.leg = scatter(REAL(ww,jj),IMAG(ww,jj),10,'o','MarkerEdgeColor','k','MarkerFaceColor','k',...
            'MarkerFaceAlpha',1,'LineWidth',1.5);

        hh(jj) = h.grand;
        
        uistack(h.grand,'top')
    end
    
%     scatter(1,0,40,'o','MarkerEdgeColor','r','MarkerFaceColor','r',...
%         'MarkerFaceAlpha',1,'LineWidth',1.5);
    
    leg = legend(hh,legLabel{ww,:});
    leg.Title.String = 'Frequency / Velocity';
    leg.Location = 'northwest';
    leg.Position = leg.Position + [-0.2 0 0 0];
    legend boxoff
end

%% Complex Gain: Normalized Bode Plot
FIG = figure (3); clf
FIG.Color = 'w';
FIG.Position = [100 100 800 700];
FIG.Name = 'Normalized Bode';
movegui(FIG,'center')
hold on

ax1 = subplot(2,1,1) ; hold on
    ax1.FontSize = 12;
    ax1.XLim = [0 12.5]; 
    ax1.YLim = [0 1.0];
    ax1.XLabel.FontSize = 14;
    ax1.XLabel.Color = 'w';
 	ax1.YLabel.String = ['Gain (' char(176) '/' char(176) ')'];
    ax1.YLabel.FontSize = ax1.XLabel.FontSize;
  	ax1.XTick = Freq;
    ax1.XTickLabel = '';

    errorbar(Freq,GAIN_NORM,2*GAIN_NORM_STD,'-b','LineWidth',4)
    
ax2 = subplot(2,1,2) ; hold on
    ax2.FontSize = ax1.FontSize;
    ax2.XLim = ax1.XLim;
    ax2.YLim = [-180 180];
   	ax2.XLabel.String = 'Frequency (Hz)';
    ax2.XLabel.FontSize = 14;
    ax2.XLabel.Color = 'k';
 	ax2.YLabel.String = ['Phase (' char(176) ')'];
    ax2.YLabel.FontSize = ax1.XLabel.FontSize;
    ax2.XTick = ax1.XTick;
    ax2.YTick = -180:60:180;
    
    errorbar(Freq,PHASE_NORM,2*PHASE_NORM_STD,'-b','LineWidth',4)
    plot([0 12.5],[0 0],'--g','LineWidth',1)
    
ax3 = axes;
    ax3.Position = ax1.Position + [0 -0.00 0 0];
    ax3.FontSize = ax1.FontSize ;
    ax3.Color = 'none';
    ax3.YAxisLocation = 'right';
    ax3.YAxis.Color = 'none';
    ax3.XAxisLocation = 'top';
    ax3.XLim = ax1.XLim;
    ax3.XTick = ax1.XTick;
    
    velLabel = cellfun(@(x) num2str(x), num2cell(vel),'UniformOutput',false);
    
    ax3.XTickLabels = velLabel;
    ax3.XLabel.String = ['Peak Velocity (' char(176) '/s)'];
    ax3.XLabel.FontSize = ax1.YLabel.FontSize;

%% Complex Gain: All on One Bode Plot
FIG = figure (4); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 4 5];
FIG.Name = 'All Bode';
movegui(FIG,'center')
hold on
clear h ax
cmap = hsv(nAmp);
for ww = 1:nAmp
    ax(1) = subplot(2,1,1) ; hold on
        ax(1).FontSize = 8;
        ax(1).XLim = [0 12.5];
        ax(1).YLim = [0 1.0];
        ax(1).XLabel.FontSize = 10;
        ax(1).XLabel.Color = 'w';
        ax(1).YLabel.String = ['Gain (' char(176) '/' char(176) ')'];
        ax(1).YLabel.FontSize = ax(1).XLabel.FontSize;
        ax(1).XTick = Freq;
        ax(1).XTickLabel = '';

        [~,h.gain(ww)] = PlotPatch(GAIN(ww,:),GAIN_STD(ww,:),Freq,1,HeadFree{ww}.N{1,1},cmap(ww,:),[0.4 0.4 0.6],0.2,2);
     	h.gain(ww).Marker = '.';
        h.gain(ww).MarkerSize = 15;

    ax(2) = subplot(2,1,2) ; hold on
        ax(2).FontSize = ax(1).FontSize;
        ax(2).XLim = ax(1).XLim;
        ax(2).YLim = [-180 180];
        ax(2).XLabel.String = 'Frequency (Hz)';
        ax(2).XLabel.FontSize = 10;
        ax(2).XLabel.Color = 'k';
        ax(2).YLabel.String = ['Phase (' char(176) ')'];
        ax(2).YLabel.FontSize = ax(1).XLabel.FontSize;
        ax(2).XTick = ax(1).XTick;
        ax(2).YTick = -180:60:180;

        [~,h.phase(ww)] = PlotPatch(PHASE(ww,:),PHASE_STD(ww,:),Freq,1,HeadFree{ww}.N{1,1},cmap(ww,:),[0.4 0.4 0.6],0.2,2);
        h.phase(ww).Marker = '.';
        h.phase(ww).MarkerSize = 15;
        plot([0 12.5],[0 0],'--k','LineWidth',1)
end
leg = legend(h.gain,string(Amp));
leg.Box = 'off';
leg.Title.String = 'Amplitude (°)';

%% Complex Gain: All Bode Plot
FIG = figure (5); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 12 5];
FIG.Name = 'All Bode';
movegui(FIG,'center')
hold on
clear h ax
cmap = hsv(nAmp);
for ww = 1:nAmp
    ax(1) = subplot(2,nAmp,ww) ; hold on
        ax(1).FontSize = 8;
        ax(1).XLim = [0 12.5];
        ax(1).YLim = [0 1.0];
        ax(1).XLabel.FontSize = 10;
        ax(1).XLabel.Color = 'w';
        ax(1).YLabel.String = ['Gain (' char(176) '/' char(176) ')'];
        ax(1).YLabel.FontSize = ax(1).XLabel.FontSize;
        ax(1).XTick = Freq;
        ax(1).XTickLabel = '';

        [~,h.gain] = PlotPatch(GAIN(ww,:),GAIN_STD(ww,:),Freq,1,1,cmap(ww,:),[0.4 0.4 0.6],0.2,2);
     	h.gain.Marker = '.';
        h.gain.MarkerSize = 15;

    ax(2) = subplot(2,nAmp,ww+nAmp) ; hold on
        ax(2).FontSize = ax(1).FontSize;
        ax(2).XLim = ax(1).XLim;
        ax(2).YLim = [-180 180];
        ax(2).XLabel.String = 'Frequency (Hz)';
        ax(2).XLabel.FontSize = 10;
        ax(2).XLabel.Color = 'k';
        ax(2).YLabel.String = ['Phase (' char(176) ')'];
        ax(2).YLabel.FontSize = ax(1).XLabel.FontSize;
        ax(2).XTick = ax(1).XTick;
        ax(2).YTick = -180:60:180;

        [~,h.phase] = PlotPatch(PHASE(ww,:),PHASE_STD(ww,:),Freq,1,1,cmap(ww,:),[0.4 0.4 0.6],0.2,2);
        h.phase.Marker = '.';
        h.phase.MarkerSize = 15;
        plot([0 12.5],[0 0],'--g','LineWidth',1)

	ax(3) = axes;
        ax(3).Position = ax(1).Position + [0 -0.00 0 0];
        ax(3).FontSize = ax(1).FontSize ;
        ax(3).Color = 'none';
        ax(3).YAxisLocation = 'right';
        ax(3).YAxis.Color = 'none';
        ax(3).XAxisLocation = 'top';
        ax(3).XLim = ax(1).XLim;
        ax(3).XTick = ax(1).XTick;

        velLabel = cellfun(@(x) num2str(x), num2cell(Vel(ww,:)),'UniformOutput',false);

        ax(3).XTickLabels = velLabel;
        ax(3).XLabel.String = ['Peak Velocity (' char(176) '/s)'];
        ax(3).XLabel.FontSize = ax(1).YLabel.FontSize;
        % set(ax,'XScale','log')
end

%% Complex Gain: All Bode Plot for EMD mean temporal frequency

emd_data = load('C:\Users\BC\Box\Research\Manuscripts\Head\EMD\data\EMD_Model=1_Wave=30_Amp=15_Accp=4.95_Delay=35.mat');

FIG = figure (5); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 14 5];
FIG.Name = 'All Bode';
movegui(FIG,'center')
hold on
clear h ax
cmap = hsv(nAmp);
mean_temp_freq = 4*Freq.*Amp / 30;
ax = gobjects(nAmp,2);
for ww = 1:nAmp
    ax(ww,1) = subplot(2,nAmp,ww) ; hold on
        ax(ww,1).YLabel.String = 'Gain (°)';
        [~,h.gain] = PlotPatch(GAIN(ww,:),GAIN_STD(ww,:), ... 
            mean_temp_freq(ww,:),1,1,cmap(ww,:),[0.4 0.4 0.6],0.2,2);
     	h.gain.Marker = '.';
        h.gain.MarkerSize = 15;
        
        plot(emd_data.mean_tempfreq_raw, emd_data.Mag_Raw, 'k', 'LineWidth', 2)
        plot(emd_data.mean_tempfreq_head, emd_data.Mag_Head, '.-b', 'LineWidth', 2, ...
                        'Marker', '.', 'MarkerSize', 15)

    ax(ww,2) = subplot(2,nAmp,ww+nAmp) ; hold on
        ax(ww,2).YLim = [-180 180];
        ax(ww,2).YLabel.String = 'Phase (°)';
        ax(ww,2).YTick = -180:60:180;

        [~,h.phase] = PlotPatch(PHASE(ww,:),PHASE_STD(ww,:), ...
            mean_temp_freq(ww,:),1,1,cmap(ww,:),[0.4 0.4 0.6],0.2,2);
        h.phase.Marker = '.';
        h.phase.MarkerSize = 15;
        plot([0 25],[0 0],'--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
        
        plot(emd_data.mean_tempfreq_raw, rad2deg(emd_data.Phase_Raw), 'k', 'LineWidth', 2)
        plot(emd_data.mean_tempfreq_head, rad2deg(emd_data.Phase_Head), '.-b', 'LineWidth', 2, ...
                        'Marker', '.', 'MarkerSize', 15)

end
set(ax, 'FontSize', 10, 'LineWidth', 1.5, 'Xlim', [0 25], 'XTick', 0:5:25)
linkaxes(ax(:,1), 'y')
linkaxes(ax(:,2), 'y')
linkaxes(ax, 'x')

XLabelHC = get(ax(:,2), 'XLabel');
set([XLabelHC{:}], 'String', 'Mean Temporal Frequency  (Hz)')

%% Gain & Phase vs Amplitude
FIG = figure (5); clf
FIG.Color = 'w';
FIG.Position = [100 100 700 700];
FIG.Name = 'Gain & Phase vs Amplitude';
movegui(FIG,'center')
hold on

cmap = jet(nFreq);
clear h
for jj = 1:nFreq
    ax = subplot(2,1,1); hold on
    ax.FontSize = 12;
	ax.YLabel.String = ['Gain (' char(176) '/' char(176) ')'];
    ax.YLabel.FontSize = 14;
    ax.XLabel.String = ['Amplitude (' char(176) ')'];
    ax.XLabel.FontSize = ax.YLabel.FontSize;
    ax.XLim = [7 19];
    ax.YLim = [0 1];
    ax.XTick = Amp;
    
    h.freq(jj) = errorbar(Amp,GAIN(:,jj),2*GAIN_STD(:,jj),'-','Color',cmap(jj,:),'LineWidth',2);
    
%     [~,h.freq(jj)] = PlotPatch(GAIN(:,jj),GAIN_STD(:,jj),Amp,1,1,cmap(jj,:),[0.5 0.5 0.6],0.1,3);

    ax = subplot(2,1,2); hold on
    ax.FontSize = 12;
	ax.YLabel.String = ['Phase (' char(176) ')'];
    ax.YLabel.FontSize = 14;
    ax.XLabel.String = ['Amplitude (' char(176) ')'];
    ax.XLabel.FontSize = ax.YLabel.FontSize;
    ax.XLim = [7 19];
    ax.YLim = [-150 150];
    ax.XTick = Amp;
    h.freq(jj) = errorbar(Amp,PHASE(:,jj),2*PHASE_STD(:,jj),'-','Color',cmap(jj,:),'LineWidth',2);
    
%     [~,h.freq(jj)] = PlotPatch(GAIN(:,jj),GAIN_STD(:,jj),Amp,1,1,cmap(jj,:),[0.5 0.5 0.6],0.1,3);
    
end

leg = legend(h.freq,num2strcell(Freq(1:jj)));
leg.Title.String = 'Frequency (Hz)';
legend boxoff

%% Complex Gain: Amplitude Change
FIG = figure (6); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 4 4];
movegui(FIG,'center')
FIG.Name = filename;
for ww = 1:nAmp
   FIG.Name = [FIG.Name '_' num2str(Amp(ww))];  
end
[ax,h] = ComplexAxes(gains,-15); hold on
ax.XLim = [-0.4 0.8];
ax.YLim = [-0.6 0.6];
% ax.XLim = [-0.04 0.1];
% ax.YLim = [-0.08 0.08];
% h.circle(end).Color = 'r';
% h.origin(1).Color = 'r';
set(h.text,'Color','k')
clear h

% for ww = 1:nAmp
%     for jj = 1:nFreq 
%         h.trial = scatter(Real{ww}(:,jj), Imag{ww}(:,jj), 15, 'o','MarkerEdgeColor','k',...
%             'MarkerFaceColor',cList(jj,:), 'MarkerFaceAlpha', 0.3, 'LineWidth', 0.5);
%     end
% end

for jj = 1:nFreq 
    h.grand(jj) = plot(REAL(:,jj),IMAG(:,jj),'-o','Color',cList(jj,:),'MarkerSize',5,'LineWidth',3.5,...
        'MarkerFaceColor','k','MarkerEdgeColor','none');
    uistack(h.grand(jj),'top')
    
%     h.text = text(1.1*REAL(:,jj),1.1*IMAG(:,jj),num2strcell(Amp));
%     set(h.text,'FontSize',8)
    
%     for ww = 1:nAmp
%         rSTD = PolarSTD(Real{ww}(:,jj),Imag{ww}(:,jj),[REAL(ww,jj) IMAG(ww,jj)]);
%         
%      	[h.std] = draw_ellipse([REAL(ww,jj) IMAG(ww,jj)], 3*rSTD, 0.5, 0, 90, cList(jj,:)); hold on
%         h.std{1}.FaceAlpha = 0.2;
%         for kk = 2:length(h.std)
%            delete(h.std{kk}) 
%         end
%         
%         scatter(REAL(ww,jj),IMAG(ww,jj),15,'o','MarkerEdgeColor','k','MarkerFaceColor',cList(jj,:),...
%             'MarkerFaceAlpha',1,'LineWidth',1.5);
%     end
end

leg = legend(h.grand,num2strcell(Freq));
leg.Title.String = 'Frequency (Hz)';
leg.Location = 'northwest';
leg.Position = leg.Position;
legend boxoff

end