% this is the script that performs the basic simulation for rev phi
% patterns and EMD response (arb. units)s

close all; clear all; % first we need to make the eye_filters and the patterns
make_eye_filters; % currently this is w. 4.6 deg ommatidia, 72 total

spatPeriod = 3.75*[8,16,24];
nPeriod = length(spatPeriod);

x_num = 96;  %% x contains rotations all around the display
for kk = 1:nPeriod
    % Pattern
    [pattern] = MakePattern_SpatFreq(spatPeriod(kk));
    temp = pattern.Pats(1,:,:);
    Pats(1,:,:,kk) = temp;
end

% make a space-time plot of the stimuli
cmap = [0 0 0; 0 1/6 0;0 2/6 0; 0 3/6 0;0 4/6 0; 0 5/6 0;0 6/6 0];
colormap(cmap);
for k = 1:nPeriod
    figure(7)
    subplot(1,3,k)
    imagesc(squeeze(Pats(1,:,:,k))');colormap(cmap);
    hold on;
    ylabel('time')
    xlabel('space');
%     if k == 1
%         title('standard motion');
%     elseif k == 4
%      title('reverse-phi motion');
%     end
end

lp_Tau_HR = 15e-3;  % time constant of the lp-filter
hp_Tau_HR = 50e-3; %%time constant of the hp filter, from Borst et al, 2003
sample_rate = 1000;

pause_time = 0.2;   num_pts_pause = pause_time*sample_rate;
OL_time = 3;      num_pts_OL = OL_time*sample_rate;

rotation_rates = [0 .5 1 2 4 8 16 32 64 96 120 192 250 300 400];
num_speeds = length(rotation_rates);
num_patterns = nPeriod;

%% simulate rotation of the visual pattern for all patterns and speeds
for j = 1:length(rotation_rates)
    
    for k = 1:num_patterns
        
        index = (j-1)*num_patterns + k          
        frame_positions(1:num_pts_pause, 1) = -1;    
        ifi = rotation_rates(j)/sample_rate; %inter-frame interv    
        frame_positions((num_pts_pause+1):(num_pts_pause+num_pts_OL),1) = ...
            mod( round( ([1:num_pts_OL] - 1)*ifi), x_num) + 1;
       % second column is constant
        frame_positions(1:(num_pts_pause + num_pts_OL),2) = k;
        tic
        [sim_data(index).eye_sample, sim_data(index).HR_Motion] = OL_arena_simulation_w_hp(eye_filt, Pats, frame_positions, sample_rate, lp_Tau_HR, hp_Tau_HR);                
        %[sim_data(index).eye_sample, sim_data(index).HR_Motion] = OL_arena_simulation(eye_filt, Pats, frame_positions, sample_rate, lp_Tau_HR);                
		toc
        
%         figure(6);
%         subplot(2,1,1); image(sim_data(index).eye_sample)
%         subplot(2,1,2); image(sim_data(index).HR_Motion)
%         
%         pause
        
        sim_data(index).HR_mean_ss = mean(sim_data(index).HR_Motion(num_pts_pause+20:end,:));
        sim_data(index).HR_mean_ts = mean(sim_data(index).HR_Motion, 2);
        sim_data(index).HR_mean_ss_avg =  mean(sim_data(index).HR_mean_ss);
        
    end
    
end

% make a log plot of mean ss response
figure(4); clf;
set(4, 'Position', [100 100 680 580],'color', 'w')
subplot(2,1,1)
temp_freq = repmat((rotation_rates*3.75), nPeriod,1 )./repmat([spatPeriod]', 1, num_speeds);
speeds = repmat((rotation_rates*3.75), nPeriod,1 );
%map = lines;
%cmap = map(linspace(2,6,nPeriod),:);

hold all
for k = 1:nPeriod
   plot(temp_freq(k,:), [sim_data(k:nPeriod:end).HR_mean_ss_avg], 'o-', 'LineWidth', 2, 'MarkerSize', 6)
end

set(gca,'xscale','log','FontSize',10,'FontName','Times');
ylim([-.58 0.58]);
set(gca, 'YTick', [-0.4 -0.2 0 0.2 0.4]);
xlim([.09 60]);
set(gca, 'XTick', [0.1, 1, 10, 50]);
set(gca, 'XTickLabel', [0, 1, 10, 50]);
set(gca,'XGrid','on','Xcolor',[.3 .3 .3])
c_axes = copyobj(gca,gcf);
set(c_axes, 'color', 'none', 'xcolor', 'k', 'xgrid', 'off', 'ycolor','k', 'ygrid','off');
xlabel('temporal frequency (Hz)','color', 'k');
ylabel('EMD response (arb. units)');
title('normal rotation');
box off

subplot(2,1,2)
% plot standard resp on linear scale
hold all
for k = 1:nPeriod
   plot(speeds(k,:), [sim_data(k:nPeriod:end).HR_mean_ss_avg], 'o-', 'LineWidth', 2, 'MarkerSize', 6)
end

set(gca,'xscale','log','FontSize',10,'FontName','Times');
set(gca,'xscale','log','FontSize',10,'FontName','Times');
ylim([-.58 0.58]);
set(gca, 'YTick', [-0.4 -0.2 0 .2 .4]);
xlim([7 2000]);
set(gca, 'XTick', [10, 100, 1000]);
set(gca, 'XTickLabel', [ 10, 100, 1000]);
set(gca,'XGrid','on','Xcolor',[.3 .3 .3])
c_axes = copyobj(gca,gcf);
set(c_axes, 'color', 'none', 'xcolor', 'k', 'xgrid', 'off', 'ycolor','k', 'ygrid','off');
xlabel('velocity (deg s^-^1)','color', 'k');
ylabel('EMD response (arb. units)','color','k');
title('normal rotation');

