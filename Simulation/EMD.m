function [EMD_ALL, Pattern_ALL, Eye_ALL] = EMD(Eye,Pattern,Pos,tt,showplot,export)
% EMD:  This simulation reconstructs the fly's visual environment and calculates
%    	optic flow based on Reichart detector
%   INPUTS:
%       Eye         :   eye object containing physiological properties
%       Pattern   	:   spatial information based on LED arena
%       Pos        	:   pattern position [pattern#]
%       tt       	:   time vector
%       showplot  	:   show plot boolean
%       export  	:   export figure boolean
%   OUTPUTS:
%       EMD_ALL    	:   
%       Pattern_ALL	:   
%       Eye_ALL    	:   
%---------------------------------------------------------------------------------------------------------------------------------
if nargin<6
    export = false; % default is off
    if nargin<5
        showplot = false; % default is off
    end
end

% Initializations for HR model
Fs  = 1/mean(diff(tt));
tau = Eye.time_constant;
h   = 1/Fs;

% Discrete low-pass filter parameters
a           = h / (tau+h);
InMat       = 5*(rand(1,Eye.n_receptor) - 0.5);
FiltMat_1   = zeros(size(InMat));

% Preallocate arrays
EMD_ALL = nan(length(Pos), Eye.n_receptor);
Pattern_ALL = nan(length(Pos), 96); % 96-pixel column
Eye_ALL = nan(length(Pos), Eye.n_receptor);
HR_Motion = nan(length(Pos), Eye.n_receptor-1);
for jj = 1:length(Pos)   
   	% CALCULATE OUTPUT FROM ELEMENTARY MOTION DETECTOR
    % (HASSENTSEIN-REICHARDT)
    
 	current_frame = Pattern.Pats(:, :, Pos(jj),1);
    % upsample by factor of 10
    for k = 1:10
        Up_frame(k:10:Eye.n_pts) = current_frame(1,1:96);
    end
    
    % Filtered eye projection
    eye_sample = Up_frame*Eye.filt;    
    
    % Compute HR motion
    InMat       = eye_sample;
    FiltMat     = a*(InMat) + (1-a)*FiltMat_1; % discrete low-pass filter
    FiltMat_1   = FiltMat;
    % HR_Motion = Va(t-tau) * Vb - Vb(t-tau) * Va
    HR_Motion(jj,:) = FiltMat(1:end-1).*InMat(2:end) - FiltMat(2:end).*InMat(1:end-1); % delay and correlate
  	
    % Prepare data for space-time plot
    EMD_ALL(jj,1:end-1) = HR_Motion(jj,:);
    Pattern_ALL(jj,:)   = current_frame(1,1:96);
    Eye_ALL(jj, :)      = eye_sample;
end

% Plot
if showplot == 1
    x = [-180:3.75:-3.75 3.75:3.75:180];
    y = [-30:3.75:-3.75 3.75:3.75:30];
    mymap = [zeros(Eye.n_receptor,1) linspace(0, 1, Eye.n_receptor)' zeros(Eye.n_receptor,1)];
    mm = 1;
    for jj = 1:length(Pos)
        pat = Pattern.Pats(:, :, Pos(jj),1);

        H = figure(50); cla; clf ; hold on
        H.Renderer = 'OpenGL';
        H.Position = [100, 100, 15*40, 15*50];

        subplot(4,1,1:2)
        imagesc(x, y, pat)
        colormap(mymap);
        xlim([-180 180])
        ylim([-30 30]);
        set(gcf, 'color', 'k');
        set(gca, 'color', 'k')
        set(gca,'xcolor','w','ycolor','w');
        set(gca,'box','off');
        freezeColors;

        subplot(4,1,3:4)
        plot(HR_Motion(jj,:),'w');
        xlim([1 Eye.n_receptor])
        ylim(max(abs(HR_Motion),[],'all')*[-1 1]);
        set(gcf, 'color', 'k');
        set(gca, 'color', 'k')
        set(gca,'xcolor','w','ycolor','w');
        xlabel('ommatidia')
        ylabel('EMD output')
        set(gca,'box','off');

        % Display time stamp
        text(30, -430, sprintf('%0.3f',tt(jj)),'fontsize',10,'color',[1 1 1])

        % Export image
        if export
            imgdirUnix = 'C:/JMM/Magno_data/7_5deg_Grnd/sim5/';
            filename = sprintf([imgdirUnix 'image%04d.jpg'], mm);
            export_fig(gcf, filename, '-q95','-nocrop');
        end
        mm = mm + 1;
    end
end

end