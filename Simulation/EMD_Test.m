
%%
clear ; close all ; clc

root.Pat = 'Q:\Box Sync\Git\Arena\Patterns';
root.Fly = 'H:\EXPERIMENTS\Experiment_Sinusoid\11.25';

% Pattern
[file.Pat, path.Pat] = uigetfile({'*.mat', 'MAT-files'}, ...
    'Select pattern file', root.Pat, 'MultiSelect','off');

load(fullfile(path.Pat,file.Pat),'pattern');

% Fly head angles (fly_2_trial_2_freq_2.mat)
[file.Fly, path.Fly] = uigetfile({'*.mat', 'MAT-files'}, ...
    'Select fly angles file', fullfile(root.Fly,'Vid','Angles'), 'MultiSelect','off');

load(fullfile(path.Fly,file.Fly),'t_v','hAngles');

load(fullfile(root.Fly,file.Fly),'t_p','data');
%%
tt = TimeVector(t_v,false);
tt = tt(1:end-1);

head.pos = interp1(t_v,hAngles,tt);
head.pos = head.pos - mean(head.pos);
stim.pos = panel2deg(data(:,2));
stim.pos = interp1(t_p,stim.pos,tt);

error.time = tt;
error.pos = stim.pos - head.pos;

figure (1) ; clf
hold on
plot(tt,stim.pos);
plot(tt,head.pos);
plot(tt,error.pos);


error.pos   = deg2panel(error.pos) + 20;
stim.pos    = deg2panel(stim.pos) + 20;

%%
% make eye filters
delta_phi = 4.6*pi/180; % angle between adjacent ommatidia
n_receptor = 64;
Eye = EYE(delta_phi,n_receptor);

[EMD_ALL, Pattern_ALL, Eye_ALL] = EMDsim_v5_7_5deg_Grnd_new(Eye,pattern,stim.pos,error.time);

[EMD_ALL_filt] = HR_sim(Eye_ALL);

EMD_ALL_filt = EMD_ALL_filt(:,1:end-1); % resize array 1-63

% make space-time plot
figure; subplot(1,2,1); imagesc(Pattern_ALL)
ylabel('time (s)'); set(gca,'YTick',[0 2.0*160])
set(gca,'YTickLabel',{'0', '2.0'})
set(gca,'XTick',[1 48 96])
set(gca,'XTickLabel',{'-180','0', '180'})
title('Intensity space-time plot');
xlabel('Pattern position'); colormap bone; freezeColors;
ylim([0 2.0*160]);

% EMD
subplot(1,2,2); 
imagesc(EMD_ALL_filt);
set(gca,'YTick',[0 2.0*160]); set(gca,'YTickLabel',{'0', '2.0'})
colormap(jet);
xlabel('ommatidia'); title('EMD'); ylim([0 2.0*160]); xlim([1 8]);
set(gca,'XTick',[1 8])
set(gca,'XTickLabel',{'1','8'})




