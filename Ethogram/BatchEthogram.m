%% Select File(s) %%
%---------------------------------------------------------------------------------------------------------------------------------
clear;close all;clc

root = 'E:\Walking_Experiments\SOS\mat';

[FILES, dirpath] = uigetfile('*.mat', 'Select fly trials', root, 'MultiSelect','on');
FILES = cellstr(FILES)';

nTrial = length(FILES); % total # of trials
playback = 5;

%% Load video data & run tracking software %%
%---------------------------------------------------------------------------------------------------------------------------------
for jj = 1:nTrial
    % Load video data
    load([dirpath FILES{jj}]); % load video data
    disp('Load File: Done')
    
    etho=MakeEthogram(rawVid,playback);
    
    

    % Save data
    disp('Save Data...')
    save(fullfile(root,'\ethograms', FILES{jj}),'etho');
    
    close all
    clc
end