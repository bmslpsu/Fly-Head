% Batch Ethogram
% Audrey Blizard
% 12/19/2019
% Loads in a series of .mat files containing fly videos. It then allows the
% user to run the MakeEthogram function on each file and save the results.

%% Select File(s) %%
%---------------------------------------------------------------------------------------------------------------------------------
clear;close all;clc

root = 'E:\Walking_Experiments\SOS\mat';                   % PC file path
%root = '/Volumes/Data_Audrey/Walking_Experiments/SOS/mat';  % Mac file path

[FILES, dirpath] = uigetfile('*.mat', 'Select Walking Videos', root, 'MultiSelect','on');
FILES = cellstr(FILES)';

nTrial = length(FILES); % total # of trials
playback = 5;           % used in function, increment of frames shown

%% Load video data & run tracking software %%
%---------------------------------------------------------------------------------------------------------------------------------
for jj = 1:nTrial
    % Load video data
    load([dirpath FILES{jj}]);              % load video data
    disp('Load File: Done')
    disp(FILES{jj})
    etho=MakeEthogram(rawVid,playback);     % calls MakeEthogram Function

    % Save data
    disp('Save Data...')
    save(fullfile(root,'\ethograms', FILES{jj}),'etho');   % PC save
    %save(fullfile(root,'/ethograms', FILES{jj}),'etho');    % Mac save
    
    close all
    clc
end

disp('All Trials Completed')