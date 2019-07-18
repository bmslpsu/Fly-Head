% this is the script that performs the basic simulation for rev phi
% patterns and EMD response (arb. units)s

clear all; % first we need to make the eye_filters and the patterns
make_eye_filters; % currently this is w. 4.6 deg ommatidia, 72 total
make_multiwidth_phi_rev_phi_pattern_grayscale; %%returns the visual pattern
lp_Tau_HR = 15e-3;  % time constant of the lp-filter
hp_Tau_HR = 50e-3; %%time constant of the hp filter, from Borst et al, 2003
sample_rate = 1000;

pause_time = .2;   num_pts_pause = pause_time*sample_rate;
OL_time = 3;      num_pts_OL = OL_time*sample_rate;

rotation_rates = [0 .5 1 2 4 8 16 32 64 96 120 192 250 300 400];
num_speeds = length(rotation_rates);
num_patterns = 6;

%% simulate rotation of the visual pattern for all patterns and speeds
for j = 1:length(rotation_rates)
    for k = 1:num_patterns
        index = (j-1)*num_patterns + k          
        frame_positions(1:num_pts_pause, 1) = -1;    
        ifi = rotation_rates(j)/sample_rate; %inter-frame interv    
        frame_positions((num_pts_pause+1):(num_pts_pause+num_pts_OL),1) = ...
            mod( round( ([1:num_pts_OL] - 1)*ifi), x_num) + 1;
       % second channel is constant
        frame_positions(1:(num_pts_pause + num_pts_OL),2) = k;
        tic
        [sim_data(index).eye_sample, sim_data(index).HR_Motion] = OL_arena_simulation_w_hp(eye_filt, Pats, frame_positions, sample_rate, lp_Tau_HR, hp_Tau_HR);        
		toc
        sim_data(index).HR_mean_ss = mean(sim_data(index).HR_Motion(num_pts_pause+20:end,:));
        sim_data(index).HR_mean_ts = mean(sim_data(index).HR_Motion, 2);
        sim_data(index).HR_mean_ss_avg =  mean(sim_data(index).HR_mean_ss);
	end
end

% make a log plot of mean ss response to phi and rev phi motion
figure(4); clf;
set(4, 'Position', [100 100 680 580],'color', 'w')
subplot(2,2,1)
temp_freq = repmat((rotation_rates*3.75), 3,1 )./repmat([90 60 30]', 1, num_speeds);
speeds = repmat((rotation_rates*3.75), 3,1 );
map = lines;
cmap = map(4:6,:);

hold all
for k = 1:3
   plot(temp_freq(k,:), [sim_data(k:6:end).HR_mean_ss_avg], 'o-', 'LineWidth', 2, 'MarkerSize', 6, 'color' ,cmap(k,:))
end

set(gca,'xscale','log','FontSize',10,'FontName','Times');
ylim([-.58 0.58]);
set(gca, 'YTick', [-0.4 -0.2 0 .2 .4]);
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

subplot(2,2,3)
% plot standard resp on linear scale
hold all
for k = 1:3
   plot(speeds(k,:), [sim_data(k:6:end).HR_mean_ss_avg], 'o-', 'LineWidth', 2, 'MarkerSize', 6, 'color',cmap(k,:))
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

% plot rev_phi on a log scale
subplot(2,2,2)
hold all
for k = 1:3
   plot(temp_freq(k,:), [sim_data((3+k):6:end).HR_mean_ss_avg], 'o-', 'LineWidth', 2, 'MarkerSize', 6, 'color',cmap(k,:))
end

set(gca,'xscale','log','FontSize',10,'FontName','Times');
ylim([-.58 0.58]);
set(gca, 'YTick', [-0.4 -0.2 0 .2 .4]);
xlim([.09 60]);
set(gca, 'XTick', [0.1, 1, 10, 50]);
set(gca, 'XTickLabel', [0, 1, 10, 50]);
set(gca,'XGrid','on','Xcolor',[.3 .3 .3])
c_axes = copyobj(gca,gcf);
set(c_axes, 'color', 'none', 'xcolor', 'k', 'xgrid', 'off', 'ycolor','k', 'ygrid','off');
xlabel('temporal frequency (Hz)','color', 'k');
ylabel('EMD response (arb. units)');
title('reverse phi rotation');
box off

subplot(2,2,4)
% plot rev_phi resp on linear scale
hold all
for k = 1:3
   plot(speeds(k,:), [sim_data((3+k):6:end).HR_mean_ss_avg], 'o-', 'LineWidth', 2, 'MarkerSize', 6, 'color',cmap(k,:))
end

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
title('reverse phi rotation');

hhh = legend('\lambda=90^o', '\lambda=60^o', '\lambda=30^o', 'Location','NorthWest');
set(hhh,'Linewidth',0.1)

for j = 1:4
subplot(2,2,j)
semilogx([.001 10000], [0 0], '--k');
end