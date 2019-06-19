function [FIG] = MakeFig_Sine_HeadFree_ComplexGain()
%% MakeFig_Sine_HeadFree_ComplexGain:
%   INPUTS:
%       -
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
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
    HeadFree{ww} = load(fullfile(root,FILES{ww}),'TRIAL','GRAND','U','N');
end
%% Complex Gain Calculations
%---------------------------------------------------------------------------------------------------------------------------------
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

cList = jet(nFreq);

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
                
                if jj==1 && I<0
                    I = abs(I);
                    R = abs(R);
                elseif jj==5 && (I>0 || R<0)
                 	I = -abs(I);
                    R = abs(R);
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

GAIN = cellfun(@(x) nanmean(x,1), Gain, 'UniformOutput', false);
PHASE = cellfun(@(x) nanmean(x,1), Phase, 'UniformOutput', false);
GAIN_STD = cellfun(@(x) nanstd(x,[],1), Gain, 'UniformOutput', false);
PHASE_STD = cellfun(@(x) nanstd(x,[],1), Phase, 'UniformOutput', false);

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

gains = 0.2:0.2:1;

%% Complex Gain: Normalized Amplitude
%---------------------------------------------------------------------------------------------------------------------------------
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
    
    h.error = plot([1 REAL_NORM(jj)],[0 IMAG_NORM(jj)],'Color','r','LineWidth',1);
    
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
%---------------------------------------------------------------------------------------------------------------------------------
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
    [ax,~] = ComplexAxes(gains);
    ax.Title.String = [num2str(Amp(ww)) char(176)];

  	for jj = 1:nFreq
        h.trial = scatter(Real{ww}(:,jj), Imag{ww}(:,jj), 40, 'o','MarkerEdgeColor','k',...
            'MarkerFaceColor',cList(jj,:), 'MarkerFaceAlpha', 0.65, 'LineWidth', 0.5);
    end
    
    for jj = 1:nFreq
        h.r = plot([0 REAL(ww,jj)],[0 IMAG(ww,jj)],'Color',[0 0 0 1],'LineWidth',1);
        
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
    
    scatter(1,0,40,'o','MarkerEdgeColor','r','MarkerFaceColor','r',...
        'MarkerFaceAlpha',1,'LineWidth',1.5);
    
    leg = legend(hh,legLabel{ww,:});
    leg.Title.String = 'Frequency / Velocity';
    leg.Location = 'northwest';
    leg.Position = leg.Position + [-0.2 0 0 0];
    legend boxoff
end

%% Complex Gain: Normalized Bode Plot
%---------------------------------------------------------------------------------------------------------------------------------
FIG = figure (3); clf
FIG.Color = 'w';
FIG.Position = [100 100 800 800];
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

%% Complex Gain: All Bode Plot
%---------------------------------------------------------------------------------------------------------------------------------
FIG = figure (4); clf
FIG.Color = 'w';
FIG.Position = [100 100 1200 800];
FIG.Name = 'All Bode';
movegui(FIG,'center')
hold on
for ww = 1:nAmp
    ax1 = subplot(2,nAmp,ww) ; hold on
        ax1.FontSize = 12;
        ax1.XLim = [0 12.5];
        ax1.YLim = [0 1.0];
        ax1.XLabel.FontSize = 14;
        ax1.XLabel.Color = 'w';
        ax1.YLabel.String = ['Gain (' char(176) '/' char(176) ')'];
        ax1.YLabel.FontSize = ax1.XLabel.FontSize;
        ax1.XTick = Freq;
        ax1.XTickLabel = '';

        errorbar(Freq,GAIN(ww,:),2*GAIN_STD(ww,:),'-b','LineWidth',4)

    ax2 = subplot(2,nAmp,ww+nAmp) ; hold on
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

        errorbar(Freq,PHASE(ww,:),2*PHASE_STD(ww,:),'-b','LineWidth',4)
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

        velLabel = cellfun(@(x) num2str(x), num2cell(Vel(ww,:)),'UniformOutput',false);

        ax3.XTickLabels = velLabel;
        ax3.XLabel.String = ['Peak Velocity (' char(176) '/s)'];
        ax3.XLabel.FontSize = ax1.YLabel.FontSize;
end

%% Gain vs Amplitude
%---------------------------------------------------------------------------------------------------------------------------------
FIG = figure (5); clf
FIG.Color = 'w';
FIG.Position = [100 100 800 800];
FIG.Name = 'Gain vs Amplitude';
movegui(FIG,'center')
hold on

cmap = jet(nFreq);

for jj = 1:nFreq
    ax = subplot(2,1,1); hold on
    ax.FontSize = 12;
	ax.YLabel.String = ['Gain (' char(176) '/' char(176) ')'];
    ax.YLabel.FontSize = 14;
    ax.XLabel.String = ['Amplitude (' char(176) ')'];
    ax.XLabel.FontSize = ax.YLabel.FontSize;
    ax.XLim = [7 19];
    ax.YLim = [0 1];
    
    h.freq(jj) = errorbar(Amp,GAIN(:,jj),GAIN_STD(:,jj),'Color',cmap(jj,:),'LineWidth',3);
    
%     [~,h.freq(jj)] = PlotPatch(GAIN(:,jj),GAIN_STD(:,jj),Amp,1,1,cmap(jj,:),[0.5 0.5 0.6],0.1,3);

    ax = subplot(2,1,2); hold on
    ax.FontSize = 12;
	ax.YLabel.String = ['Phase (' char(176) ')'];
    ax.YLabel.FontSize = 14;
    ax.XLabel.String = ['Amplitude (' char(176) ')'];
    ax.XLabel.FontSize = ax.YLabel.FontSize;
    ax.XLim = [7 19];
    ax.YLim = [-150 150];
    
    h.freq(jj) = errorbar(Amp,PHASE(:,jj),PHASE_STD(:,jj),'Color',cmap(jj,:),'LineWidth',3);
    
%     [~,h.freq(jj)] = PlotPatch(GAIN(:,jj),GAIN_STD(:,jj),Amp,1,1,cmap(jj,:),[0.5 0.5 0.6],0.1,3);
    
end

leg = legend(h.freq,num2strcell(Freq));
leg.Title.String = 'Frequency (Hz)';
legend boxoff

end