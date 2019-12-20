% Plot Ethogram
% Audrey Blizard
% 12/19/2019
% Loads in a series of ethogram vectors into large matrix. Outputs image
% with greyscale ethogram and colored bars for each fly.

%% Select File(s)
clear;close all;clc

root = 'S:\Public\Audrey\Walking SOS\Raw Data\mat\ethograms';   % PC file path
%root = '/Volumes/Data_Audrey/Walking_Experiments/SOS/mat';     % Mac file path

[D,I,N,U,T,FILES,PATH,basename] = GetFileData(root);
nFrame = 1680;                                                  % frames in 20s of data
percent = zeros(N{1,3}, 1);
ethofull = NaN(N{1,3}, nFrame);


%% Load ethogram data & create matrix 

for jj = 1:N{1,3}
    load([PATH FILES{jj}]);                                 % loads individual ethogram
    ethoshort=etho(1:nFrame)';                              % cuts down and tranposes ethogram
    percent(jj,1) = nnz(ethoshort)/length(ethoshort);       % finds percent of time walking
    ethofull(jj,:) = ethoshort;                             % adds ethogram to matrix
end

%% Process Matrix

ethofly = cat(2,ethofull+1, repmat(I{:,1}+2,1,.1*nFrame));  % adds fly columns to matrix and scales values
[~,idxpercent] = sort(percent);                             % index to sort by percent
ethosort = ethofull(idxpercent,:);                          % sorts ethogram by percent
ethoflysort = ethofly(idxpercent,:);                        % sorts ethogram with flies by percent
cmap = colormap(jet(22));                                   % creates fly colormap
cmap = cat(1,[.8 .8 .8], [0 0 0],  cmap);                   % adds colors for ethogram to colormap
sizedflysort = repelem(ethoflysort, 11,1);              % repeates each row 11 times for sizing
imshow(sizedflysort, cmap);                             % displays chart as image

%% Create Box Plot

percentbyfly = cell(1,N{1,1});
for kk = 1:N{1,3}
    flyidx = I{kk,1};
    percentbyfly{1,flyidx} = [percentbyfly{1,flyidx};percent(kk,1)];
end

fly.mean = cellfun(@(x) mean(x), percentbyfly, 'uniformoutput', true);
fly.median = cellfun(@(x) median(x), percentbyfly, 'uniformoutput', true);
fly.std = cellfun(@(x) std(x), percentbyfly, 'uniformoutput', true);
grand.mean = mean(fly.mean);
grand.median = median(fly.median);
grand.stdmean = std(fly.mean);
grand.stdmedian = std(fly.median);

FIG = figure (2); clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [1 1 3 3];
movegui(FIG,'center')
ax = gca;
hold on
boxplot(fly.mean.*100); 

FontSize = 8;
ylabel('Percent Walking');
ax.YLim = [0,100]




