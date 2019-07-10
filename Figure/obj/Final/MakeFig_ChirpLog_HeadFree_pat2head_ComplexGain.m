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
catIdx = 5;
xIdx = 1;

fIdx        = 11:200;
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


GAIN_STD  = cellfun(@(x) std(x,[],2), Gain, 'UniformOutput', false);
GAIN_STD  = cat(2,GAIN_STD{:});
PHASE_STD = cellfun(@(x) circ_std(x,[],[],2), Phase, 'UniformOutput', false);
PHASE_STD = cat(2,PHASE_STD{:});

% GAIN_r  = abs(REAL + 1i*IMAG);
% PHASE_r = angle(REAL + 1i*IMAG);

velIdx = fliplr([1 23 44 64 length(fIdx)]);

cmap = jet(length(fIdx));

gains = 0.1:0.1:0.2;
gains = 0.2:0.2:1;

%% Complex Gain: Normalized amplitudes
%---------------------------------------------------------------------------------------------------------------------------------
FIG = figure (1); clf
FIG.Color = 'w';
FIG.Position = [100 100 700 700];
FIG.Name = filename;
movegui(FIG,'center')
hold on

[ax,~] = ComplexAxes(gains);
ax.Colormap = cmap;

% for jj = 1
%     for kk = 1:size(CmplxGain{jj},1)
%         R = Real{jj}(kk,:);
%         I = Imag{jj}(kk,:);
%         h.trial = scatter(R, I, 10, cmap(kk,:),'o','MarkerEdgeColor', cmap(kk,:),...
%             'MarkerFaceColor', cmap(kk,:), 'MarkerFaceAlpha', 0.5, 'LineWidth', 0.1);
%     end
% end

for jj = 1:nAmp
    for kk = 1:size(CmplxGain{jj},1)
        if (kk>=velIdx(jj+1)) && (kk<=velIdx(jj))
            h.grand = scatter(REAL(kk,jj), IMAG(kk,jj), 30, cmap(kk,:), 'o', 'MarkerEdgeColor', 'k',...
                'MarkerFaceColor', cmap(kk,:), 'MarkerFaceAlpha', 1, 'LineWidth', 0.3);
        end
    end
end
cbar = colorbar;
cbar.FontSize = 14;
cLabel = [0.3, 2:2:10];
cbar.Ticks = maptorange(cLabel,fRange,[0 1]);
cbar.TickLabels = num2strcell(cLabel);
cbar.Label.String = 'Frequency (Hz)';
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

%% Complex Gain: one amplitude
%---------------------------------------------------------------------------------------------------------------------------------
FIG = figure (10); clf
FIG.Color = 'w';
FIG.Position = [100 100 700 700];
FIG.Name = filename;
movegui(FIG,'center')
hold on

amp = 3;

[ax,~] = ComplexAxes(gains);
ax.Colormap = cmap;
ax.Title.String = [num2str(HeadFree.U{1,3}{1}(amp)) char(176)];

for kk = 1:size(CmplxGain{amp},1)
    h.grand = scatter(REAL(kk,amp), IMAG(kk,amp), 30, cmap(kk,:), 'o', 'MarkerEdgeColor', 'k',...
        'MarkerFaceColor', cmap(kk,:), 'MarkerFaceAlpha', 1, 'LineWidth', 0.3);
end

cbar = colorbar;
cbar.FontSize = 14;
cLabel = [0.3, 2:2:10];
cbar.Ticks = maptorange(cLabel,fRange,[0 1]);
cbar.TickLabels = num2strcell(cLabel);
cbar.Label.String = 'Frequency (Hz)';
% cbar.Location = 'westoutside';
% cbar.Position = cbar.Position + [-0.1 -0.2 0 0];

%% BODE
%---------------------------------------------------------------------------------------------------------------------------------
FIG = figure (3); clf
FIG.Color = 'w';
FIG.Position = [100 100 1200 700];
FIG.Name = filename;
movegui(FIG,'center')
hold on
clear ax

for jj = 1:nAmp
    ax.Gain = subplot(2,nAmp,jj); hold on
    ax.Gain.FontSize = 12;
    ax.Gain.YLabel.String = ['Gain (' char(176) '/' char(176) ')'];
    ax.Gain.YLabel.FontSize = 14;
    ax.Gain.YLim = [0 1.2];
    ax.Gain.XLim = fRange;
    ax.Gain.XTick = unique(sort([min(ax.Gain.XLim) ax.Gain.XTick]));
   	ax.Gain.XTickLabels = '';
    vel = round(Amp(jj)*2*pi*ax.Gain.XTick);
    velLabel = num2strcell(vel);

    PlotPatch(GAIN(:,jj),GAIN_STD(:,jj),Freq,3,HeadFree.N{1,1},[0 0 0.7],[0.5 0.5 0.5],0.5,2);

    ax.Phase = subplot(2,nAmp,jj+nAmp); hold on
    ax.Phase.FontSize = 12;
    ax.Phase.YLabel.String = ['Phase (' char(176) ')'];
    ax.Phase.YLabel.FontSize = 14;
    ax.Phase.XLabel.String = 'Frequency (Hz)';
    ax.Phase.XLim = fRange;
    ax.Phase.XTick = ax.Gain.XTick;
    
    PlotPatch(PHASE(:,jj),PHASE_STD(:,jj),Freq,3,HeadFree.N{1,1},[0 0 0.7],[0.5 0.5 0.5],0.5,2);
    
    ax.Vel = axes;
    ax.Vel.Position = ax.Gain.Position + [0 -0.00 0 0];
    ax.Vel.FontSize = ax.Gain.FontSize ;
    ax.Vel.Color = 'none';
    ax.Vel.YAxisLocation = 'right';
    ax.Vel.YAxis.Color = 'none';
    ax.Vel.XAxisLocation = 'top';
    ax.Vel.XLim = ax.Gain.XLim;
    ax.Vel.XTick = ax.Gain.XTick;
    ax.Vel.XTickLabels = velLabel;
    ax.Vel.XLabel.String = ['Peak Velocity (' char(176) '/s)'];
    ax.Vel.XLabel.FontSize = ax.Gain.YLabel.FontSize;
    
end

end