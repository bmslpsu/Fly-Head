function [EMD_ALL, Pattern_ALL, Eye_ALL] = EMDsim_v5_7_5deg_Grnd_new(Eye,Pattern,Pos,tt)
% This simulation reconstructs the fly's visual environment and calculates
% optic flow based on Reichart detector
showplot = false;
export = false;

% initializations for HR model
samp_rate = 1/mean(diff(tt));
tau = 0.035;
h = 1/samp_rate;

% discrete low-pass filter parameters
a = h / (tau+h);

InMat = 5*(rand(1,Eye.n_receptor) - 0.5);
FiltMat_1 = zeros(size(InMat));

% preallocate array
EMD_ALL = nan(length(Pos), Eye.n_receptor);
Pattern_ALL = nan(length(Pos), 96); % 96-pixel column
Eye_ALL = nan(length(Pos), Eye.n_receptor);

mm = 1;
for jj = 1:length(Pos)
 	pat = Pattern.Pats(:, :, Pos(jj),4);
    
   	% CALCULATE OUTPUT FROM ELEMENTARY MOTION DETECTOR
    % (HASSENTSEIN-REICHARDT)
    
 	current_frame = (Pattern.Pats(:, :, Pos(jj),8));
    % upsample by factor of 10
    for k = 1:10
        Up_frame(k:10:Eye.n_pts) = current_frame(1,1:96);
    end
    
    % get filtered eye projection
    eye_sample = Up_frame*Eye.filt;    
    
    % compute HR motion
    InMat = eye_sample;
    FiltMat = a*(InMat) + (1-a)*FiltMat_1; % discrete low-pass filter
    FiltMat_1 = FiltMat;
    % HR_Motion = Va(t-tau) * Vb - Vb(t-tau) * Va
    HR_Motion = (FiltMat(1:end-1).*InMat(2:end) - FiltMat(2:end).*InMat(1:end-1)); % delay and correlate
  	
    % prepare data for space-time plot
    EMD_ALL(jj,1:end-1) = HR_Motion;
    Pattern_ALL(jj,:) = current_frame(1,1:96);
    Eye_ALL(jj, :) = eye_sample;
    
	if showplot == 1
        H = figure(50); cla; clf
        set(H, 'Renderer','OpenGL');
        set(H, 'Position',[100, 100, 16*40, 16*50]);

        x = [-180:3.75:-3.75 3.75:3.75:180];
        y = [-30:3.75:-3.75 3.75:3.75:30];
        hold on;
        mymap = [zeros(64,1) linspace(0, 1, 64)' zeros(64,1)]; % custom color map

        subplot(4,1,1:2)
        imagesc(x, y, pat) % show display
%         subimage(x, y, pat, mymap) % show display
        colormap(mymap);
        xlim([-180 180]); ylim([-30 30]);
        set(gca,'XTick',[-180 0 180]); set(gca,'YTick',[-30 0 30]);
        set(gca, 'YTickLabel', {'30','0', '-30'})
        set(gcf, 'color', 'k');
        set(gca, 'color', 'k')
        set(gca,'xcolor','w','ycolor','w');
        set(gca,'box','off');
        freezeColors;

        subplot(4,1,3:4)
        plot(HR_Motion,'w');
        xlim([1 64]); ylim([-300 300]);
        set(gca,'XTick',[1 64]); set(gca,'YTick',[-300 0 300]);
        set(gca, 'YTickLabel', {'','', ''})
        set(gca, 'XTickLabel', {'1','64'})
        set(gcf, 'color', 'k');
        set(gca, 'color', 'k')
        set(gca,'xcolor','w','ycolor','w');
        xlabel('ommatidia')
        ylabel('EMD output')
        set(gca,'box','off');
        % display time stamp
        text(30, -430, sprintf('%0.3f',tt(jj)),'fontsize',10,'color',[1 1 1])
    end

    % export image
    if export == 1
        imgdirUnix = ['C:/JMM/Magno_data/7_5deg_Grnd/sim5/'];
        filename = sprintf([imgdirUnix 'image%04d.jpg'], mm);
        export_fig(gcf, filename, '-q95','-nocrop');
    end
    mm = mm + 1;
end

end