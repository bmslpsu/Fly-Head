% Plot Ethogram
% Audrey Blizard
% 12/19/2019
% Loads in a series of ethogram vectors. 

%% Select File(s) %%
%---------------------------------------------------------------------------------------------------------------------------------
clear;close all;clc

root = 'E:\Walking_Experiments\SOS\mat\ethograms';                   % PC file path
%root = '/Volumes/Data_Audrey/Walking_Experiments/SOS/mat';  % Mac file path

[D,I,N,U,T,FILES,PATH,basename] = GetFileData(root);
nFrame = 1680;
percent = zeros(N{1,3}, 1);
ethofull = NaN(N{1,3}, nFrame);


%% Load video data & run tracking software %%
%---------------------------------------------------------------------------------------------------------------------------------
for jj = 1:N{1,3}
    load([PATH FILES{jj}]);              % load video data
    ethoshort=etho(1:nFrame)';
    percent(jj,1) = nnz(ethoshort)/length(ethoshort);
    ethofull(jj,:) = ethoshort;
    ethocolor(jj,:) = ethoshort*I{jj,1};
end

%% Process
imagesc(ethofull)

ethofly = cat(2,22*ethofull, repmat(2*I{:,1},1,.1*nFrame));
[~,idxpercent] = sort(percent);
ethosort = ethofull(idxpercent,:);
ethocolorsort=ethocolor(idxpercent,:);
ethoflysort = ethofly(idxpercent,:);
