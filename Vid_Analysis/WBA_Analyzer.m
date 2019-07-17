%% Select File(s) %%
%---------------------------------------------------------------------------------------------------------------------------------
clear;close all;clc

root = 'H:\EXPERIMENTS\Experiment_SOS_v2\';
[FILES, PATH] = uigetfile({'*.mat', 'DAQ-files'},'Select files', root, 'MultiSelect','on');
FILES = cellstr(FILES)';

nTrial = length(FILES); % total # of trials

%% Load video data & run tracking software %%
%---------------------------------------------------------------------------------------------------------------------------------

load([PATH FILES{1}],'vidData')

vid = squeeze(vidData);
dim = size(vid);
med = median(vid,3);
%%
FIG = figure (1) ; clf
FIG.Color = [0.2 0.2 0.2];
FIG.Units = 'inches';
ax = gca;
imshow(med)
axis([0 dim(2) 0 dim(1)])
axis equal
FIG.Position = [3 3 6 4];

[Mask] = MakeWingMask(med);

%%
clc
[Wing] = WingTracker_Edge(med,Mask);
