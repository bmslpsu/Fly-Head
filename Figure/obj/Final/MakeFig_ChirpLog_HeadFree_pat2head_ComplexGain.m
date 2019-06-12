function [FIG] = MakeFig_ChirpLog_HeadFree_pat2head_ComplexGain()
%% MakeFig_ChirpLog_HeadFree_pat2head_ComplexGain:
%   INPUTS:
%       -
%   OUTPUTS:
%       FIG     : figure handle
%---------------------------------------------------------------------------------------------------------------------------------
root = 'H:\DATA\Rigid_Data\';

[CHIRP,~] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select chirp file', root, 'MultiSelect','on');
CHIRP = cellstr(CHIRP)';

HeadFree = load(fullfile(root,CHIRP{1}),'TRIAL','U','N');
%%
clearvars -except HeadFree
filename = 'ChirpLog_HeadFree_pat2head_ComplexGain';
catIdx = 6;
xIdx = 1;

fIdx      = 7:200;
Amp         = HeadFree.U{1,3}{1};
nAmp        = HeadFree.N{1,3};
Freq        = HeadFree.TRIAL{1,1}{1,catIdx}.BodeFv(fIdx,xIdx);
fRange      = round([min(Freq) max(Freq)],1);
Vel         = (Amp*2*pi*Freq')';
vRange      = nan(2,nAmp);
CmplxGain   = cell(nAmp,1);
for jj = 1:nAmp
    pp = 1;
    for kk = 1:HeadFree.N{1,1}
        for ii = 1:size(HeadFree.TRIAL{kk,jj},1)
            CmplxGain{jj}(:,pp) = HeadFree.TRIAL{kk,jj}{ii,catIdx}.CmplxGain(fIdx,xIdx);
            pp = pp + 1;
        end
    end
    vRange(:,jj) = [min(Vel(:,jj)) max(Vel(:,jj))];
end

Real = cellfun(@(x) real(x), CmplxGain, 'UniformOutput', false);
Imag = cellfun(@(x) imag(x), CmplxGain, 'UniformOutput', false);

REAL = cellfun(@(x) mean(x,2), Real, 'UniformOutput', false);
REAL = cat(2,REAL{:});
IMAG = cellfun(@(x) mean(x,2), Imag, 'UniformOutput', false);
IMAG = cat(2,IMAG{:});

Gain  = cellfun(@(x) abs(x), CmplxGain, 'UniformOutput', false);
Phase = cellfun(@(x) angle(x), CmplxGain, 'UniformOutput', false);

GAIN  = cellfun(@(x) mean(x,2), Gain, 'UniformOutput', false);
GAIN  = cat(2,GAIN{:});
PHASE = cellfun(@(x) circ_mean(x,[],2), Phase, 'UniformOutput', false);
PHASE = cat(2,PHASE{:});

GAIN_r  = abs(REAL + 1i*IMAG);
PHASE_r = angle(REAL + 1i*IMAG);

velIdx = [23 44 64 length(fIdx)];

cmap = jet(length(fIdx));

%% Complex Gain: One amplitudes
%---------------------------------------------------------------------------------------------------------------------------------
FIG = figure (1); clf
FIG.Color = 'w';
FIG.Position = [100 100 700 700];
FIG.Name = filename;
movegui(FIG,'center')
hold on

gains = 0.2:0.2:1;
[ax,~] = ComplexAxes(gains);
ax.Colormap = cmap;

% for jj = 2
%     for kk = 1:size(CmplxGain{jj},1)
%         R = Real{jj}(kk,:);
%         I = Imag{jj}(kk,:);
%         h.trial = scatter(R, I, 10, cmap(kk,:),'o','MarkerEdgeColor', cmap(kk,:),...
%             'MarkerFaceColor', cmap(kk,:), 'MarkerFaceAlpha', 0.5, 'LineWidth', 0.1);
%     end
% end

for jj = 1
    for kk = 1:size(CmplxGain{jj},1)
        if kk<=velIdx(end)
            h.grand = scatter(REAL(kk,jj), IMAG(kk,jj), 30, cmap(kk,:), 'o', 'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', cmap(kk,:), 'MarkerFaceAlpha', 1, 'LineWidth', 0.3);
        end
    end
end
cbar = colorbar;
cbar.FontSize = 14;
cbar.Ticks = [0 1];
cbar.TickLabels = num2strcell(round([min(Freq) max(Freq)],1));
cbar.Label.String = 'Frequency (Hz)';
% cbar2 = cbar;
% cbar2.AxisLocation = 'in';
% cbar.Location = 'westoutside';

%% Complex Gain: All amplitudes
%---------------------------------------------------------------------------------------------------------------------------------
FIG = figure (2); clf
FIG.Color = 'w';
FIG.Position = [100 100 1200 700];
FIG.Name = filename;
movegui(FIG,'center')
hold on

for jj = 1:nAmp
    subplot(2,2,jj); hold on
    [ax,~] = ComplexAxes(gains);
    ax.Colormap = cmap;
    ax.Title.String = [num2str(HeadFree.U{1,3}{1}(jj)) char(176)];

    for kk = 1:size(CmplxGain{jj},1)
        h.grand = scatter(REAL(kk,jj), IMAG(kk,jj), 30, cmap(kk,:), 'o', 'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', cmap(kk,:), 'MarkerFaceAlpha', 1, 'LineWidth', 0.3);
    end
    
    if jj==1
        cbar = colorbar;
        cbar.FontSize = 14;
        cLabel = [0.3, 2:2:10];
        cbar.Ticks = maptorange(cLabel,fRange,[0 1]);
        cbar.TickLabels = num2strcell(cLabel);
        cbar.Label.String = 'Frequency (Hz)';
        cbar.Location = 'westoutside';
        cbar.Position = cbar.Position + [-0.1 -0.2 0 0];       
    end
end

%% BODE
%---------------------------------------------------------------------------------------------------------------------------------
FIG = figure (3); clf
FIG.Color = 'w';
FIG.Position = [100 100 700 700];
FIG.Name = filename;
movegui(FIG,'center')
hold on
clear ax

ax.Gain = subplot(2,1,1); hold on
ax.Gain.FontSize = 12;
ax.Gain.YLabel.String = ['Gain (' char(176) '/' char(176) ')'];
ax.Gain.YLabel.FontSize = 14;
ax.Gain.XLim = fRange;

plot(Freq,GAIN(:,3))

plot(Freq,GAIN_r(:,3))


ax.Phase = subplot(2,1,2); hold on
ax.Phase.FontSize = 12;
ax.Phase.YLabel.String = ['Phase (' char(176) ')'];
ax.Phase.YLabel.FontSize = 14;
ax.Phase.XLabel.String = 'Frequency (Hz)';
ax.Phase.XLim = fRange;
% [b,a] = butter(2,0.5,'low');
plot(Freq,PHASE(:,3))
plot(Freq,PHASE_r(:,3))




end