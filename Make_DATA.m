function [PAT,WING,HEAD,FD,T,n,unq] = Make_DATA(rootdir)
%% Make_DATA: Reads in all raw trials, transforms data, and saves in organized structure for use with figure functions
%   INPUTS:
%       root    : root directory
%   OUTPUTS:
%       PAT     : pattern structure
%       WING   	: wings structure
%       HEAD  	: head structure
%       FD      : file data
%       T       : fly data table
%       n       : field #'s 
%       unq   	: unique fields
%---------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE INPUT %
% rootdir = 'E:\Experiment_HeadExcitation\Chirp\HeadFree\';
%---------------------------------------------------------------------------------------------------------------------------------
%% Setup Directories %%
%---------------------------------------------------------------------------------------------------------------------------------
root.daq = rootdir;
root.ang = [root.daq '\Vid\Angles\'];

% Select files
[FILES, PATH.ang] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root.ang, 'MultiSelect','on');
FILES = FILES';

PATH.daq = root.daq;

%% Process File Data %%
%---------------------------------------------------------------------------------------------------------------------------------
% Read in data from head file names
n.Trial     = length(FILES);        % total # of trials
FD.Fly      = zeros(n.Trial,1); 	% fly #
FD.Trial    = zeros(n.Trial,1);     % trial #
FD.Amp      = zeros(n.Trial,1);     % amplitude [deg]
for jj = 1:n.Trial
    temp = textscan(char(FILES{jj}), '%s', 'delimiter', '_.'); temp = temp{1} ; % read individual strings into temp variable
    FD.Fly(jj,1)    = str2double(temp{2}); % store fly #
    FD.Trial(jj,1)  = str2double(temp{4}); % store trial #
    
    if ~isnan(str2double(temp{7}))
        FD.Amp(jj,1) = str2double([temp{6} '.' temp{7}]); % store amplitude
    else
        FD.Amp(jj,1) = str2double(temp{6});
    end
end
clear temp jj
%% Set up indexing convention for files %%
% Normalize to start at fly#1 and increment by 1 for each fly & trial
%---------------------------------------------------------------------------------------------------------------------------------
unq.Fly = sort(unique(FD.Fly)); % original # of each unique fly
n.Fly = length(unq.Fly); % # flies
FD.trialFly = cell(n.Fly,2);
for kk = 1:n.Fly
    FD.trialFly{kk,1} = unq.Fly(kk);  % fly # label
	FD.trialFly{kk,2} = length(find(FD.Fly==unq.Fly(kk))); % # trials per fly
end
% Make indexing array for flies
FD.newFly = (1:n.Fly)';
FD.idxFly = zeros(n.Fly,1);
pp = 1;
for kk = 1:n.Fly
    FD.idxFly(pp:pp+FD.trialFly{kk,2}-1,1) = FD.newFly(kk); % fly index
    pp = pp + FD.trialFly{kk,2};
end
% Make indexing array for trials
pp = 1;
FD.idxTrial = zeros(n.Trial,1);
for kk = 1:n.Fly
   FD.idxTrial(pp:pp+FD.trialFly{kk,2}-1) = 1:FD.trialFly{kk,2}; % trial index
   pp = pp+FD.trialFly{kk,2};
end
% Make indexing array for amplitudes
unq.Amp 	= sort(unique(FD.Amp));	% find all unique amplitudes
n.Amp       = length(unq.Amp);  	% # of unique amplitudes
FD.idxAmp   = zeros(n.Trial,1);     % amplitude index
for kk = 1:n.Amp
   FD.idxAmp(FD.Amp == unq.Amp(kk)) = kk; % frequency index
end

fprintf('Total Flies''s: %i \n',n.Fly) ; fprintf('Total Trials''s: %i \n',n.Trial)
T = cell2table(FD.trialFly,'VariableNames',{'Fly','Trials'});
disp(T)
clear pp kk
%% Get Data %%
%---------------------------------------------------------------------------------------------------------------------------------
disp('Loading Data...')
% Preallocate data cells
for kk = 1:n.Fly
    for jj = 1:n.Amp
        % PATTERN Data
        PAT.Pos{kk,1}{jj,1}         = [];
     	PAT.Time{kk,1}{jj,1}        = [];
        PAT.Vel{kk,1}{jj,1}         = [];
        PAT.VelMean{kk,1}{jj,1} 	= [];
        PAT.VelSTD{kk,1}{jj,1}      = [];
        PAT.Freq{kk,1}{jj,1}        = [];
        PAT.Mag{kk,1}{jj,1}     	= [];
        PAT.Phase{kk,1}{jj,1}       = [];
        % WING Data
        WING.Time{kk,1}{jj,1}       = [];
        WING.Pos{kk,1}{jj,1}        = [];
        WING.Vel{kk,1}{jj,1}        = [];
        WING.VelMean{kk,1}{jj,1} 	= [];
        WING.VelSTD{kk,1}{jj,1} 	= [];
        WING.Freq{kk,1}{jj,1}       = [];
        WING.Mag{kk,1}{jj,1}     	= [];
        WING.Phase{kk,1}{jj,1}      = [];
        WING.COHR.Freq{kk,1}{jj,1}  = [];
        WING.COHR.Mag{kk,1}{jj,1}   = [];
        WING.GAIN{kk,1}{jj,1}       = [];
        WING.PHASE{kk,1}{jj,1}      = [];
        % HEAD Data
        HEAD.Time{kk,1}{jj,1}       = [];
        HEAD.Pos{kk,1}{jj,1}        = [];
        HEAD.Vel{kk,1}{jj,1}        = [];
        HEAD.VelMean{kk,1}{jj,1}    = [];
        HEAD.VelSTD{kk,1}{jj,1}     = [];
        HEAD.Freq{kk,1}{jj,1}       = [];
        HEAD.Mag{kk,1}{jj,1}      	= [];
        HEAD.Phase{kk,1}{jj,1}      = [];
    	HEAD.Err.Pos{kk,1}{jj,1}	= [];
        HEAD.Err.Freq{kk,1}{jj,1}	= [];
        HEAD.Err.Mag{kk,1}{jj,1}	= [];
        HEAD.Err.Phase{kk,1}{jj,1}	= [];
    	HEAD.COHR.Freq{kk,1}{jj,1}  = [];
        HEAD.COHR.Mag{kk,1}{jj,1}   = [];
        HEAD.GAIN{kk,1}{jj,1}       = [];
        HEAD.PHASE{kk,1}{jj,1}      = [];
    end
end
% Store data in organized cells
for kk = 1:n.Trial
    disp(kk)
    % Load HEAD & DAQ data %
    data = [];
	load([PATH.daq   FILES{kk}],'data','t_p'); % load pattern x-position
    load([PATH.ang   FILES{kk}],'hAngles','t_v'); % load head angles % time arrays
    %-----------------------------------------------------------------------------------------------------------------------------
    % Check WBF
	wing.f = 100*(data(:,6)); % wing beat frequency
    if min(wing.f)<150 || mean(wing.f)<180 % check WBF, if too low then don't use trial
        fprintf('Low WBF: Fly %i Trial %i \n',FD.Fly(kk),FD.Trial(kk))
        continue
    end
	%-----------------------------------------------------------------------------------------------------------------------------
	% Get pattern data from DAQ %
    pat.Time        = t_p; % pattern time from DAQ [s]
    pat.n           = length(pat.Time); % # of samples for pat data
    pat.Fs          = 1/mean(diff(pat.Time)); % pattern sampling frequency [Hz]
    pat.Pos         = panel2deg(data(:,2));  % pattern x-pos: subtract mean and convert to deg [deg]
	pat.Vel         = [diff(pat.Pos)/(1/pat.Fs) ; 0]; % pattern velocity [deg/s]
    pat.VelMean     = mean(abs(pat.Vel)); % mean pattern velocity [deg/s]
	pat.VelSTD      = mean(abs(pat.Vel)); % STD pattern velocity [deg/s]
    %-----------------------------------------------------------------------------------------------------------------------------
    % Get head data %
	head.Time       = t_v; % store head time vector [s]
    head.n          = length(head.Time); % # of samples for wing data
    head.Fs         = 1/mean(diff(head.Time)); % sampling frequency [Hz]
    head.Fc         = 20; % cutoff frequency [Hz]
    [b,a]           = butter(2,head.Fc/(head.Fs/2),'low'); % 2nd-order low-pass butterworth filter
    head.Pos        = filtfilt(b,a,hAngles); % filter head position [deg]
    head.Pos        = head.Pos - mean(head.Pos); % subtract DC component [deg]
    head.Vel        = filtfilt(b,a,[diff(head.Pos)./diff(head.Time) ; 0]); % calculte hea vecloity and filter again [deg/s]
    head.VelMean    = mean(abs(head.Vel)); % mean head velocity [deg/s]
    head.VelSTD     = std(abs(head.Vel)); % STD pattern velocity [deg/s]
  	%-----------------------------------------------------------------------------------------------------------------------------
    % Get wing data from DAQ %
    wing.Time       = t_p; % wing time [s]
    wing.n          = length(wing.Time); % # of samples for wing data
    wing.Fs         = 1/mean(diff(wing.Time)); % sampling frequency [Hz]
    wing.Fc         = 20; % cutoff frequency [Hz]
    [b,a]           = butter(2,wing.Fc/(wing.Fs/2)); % butterworth filter
	wing.Left       = filtfilt(b,a,(data(:,4))); % left wing [V]
    wing.Right      = filtfilt(b,a,(data(:,5))); % right wing [V]
    wing.Pos        = wing.Left - wing.Right; % dWBA (L-R) [V]
  	wing.Pos        = wing.Pos - mean(wing.Pos); % subtract mean [V]
  	wing.Vel        = [diff(wing.Pos)./wing.Fs ; 0]; % dWBA velocity [V/s]
    wing.VelMean    = mean(abs(wing.Vel)); % mean dWBA velocity [V/s]
    wing.VelSTD     = std(abs(wing.Vel)); % STD dWBA velocity [V/s]
    %-----------------------------------------------------------------------------------------------------------------------------
	% Decimate DAQ data to match head Fs %
	pat.Time        = resample(pat.Time,head.n,pat.n);
    pat.Pos         = resample(pat.Pos,head.n,pat.n);
  	pat.Vel         = resample(pat.Vel,head.n,pat.n);
    wing.Time       = resample(wing.Time,head.n,wing.n);
    wing.Pos        = resample(wing.Pos,head.n,wing.n);
    wing.Left       = resample(wing.Left,head.n,wing.n);
	wing.Right      = resample(wing.Left,head.n,wing.n);
  	wing.Vel        = resample(wing.Vel,head.n,wing.n);
    wing.f          = resample(wing.f,head.n,wing.n);
	%-----------------------------------------------------------------------------------------------------------------------------
 	head.Err.Pos    = pat.Pos - head.Pos; % calculate Error between head & pattern (retinal slip) [deg]
    %-----------------------------------------------------------------------------------------------------------------------------
    % Convert head, wings, & pattern data to frequency domain using FFT %
    [head.Freq , head.Mag , head.Phase]             = FFT(head.Time,head.Pos);
    [wing.Freq , wing.Mag , wing.Phase]             = FFT(wing.Time,wing.Pos);
    [pat.Freq  , pat.Mag  , pat.Phase ]             = FFT(pat.Time,pat.Pos);
	[head.Err.Freq , head.Err.Mag, head.Err.Phase]  = FFT(head.Time,head.Err.Pos);
    %-----------------------------------------------------------------------------------------------------------------------------
    % Calculate coherence %
    [head.cohr.mag,head.cohr.f] = mscohere(pat.Pos , head.Pos ,[],[] , head.Freq , head.Fs);
	[wing.cohr.mag,wing.cohr.f] = mscohere(pat.Pos , wing.Pos ,[],[] , wing.Freq , wing.Fs);
    %-----------------------------------------------------------------------------------------------------------------------------
	% Calculate BODE gain & phase difference for head & wings %
    head.GAIN   = head.Mag./pat.Mag;
    head.PHASE  = pat.Phase - head.Phase;
 	wing.GAIN   = wing.Mag./pat.Mag;
    wing.PHASE  = pat.Phase - wing.Phase;
    %-----------------------------------------------------------------------------------------------------------------------------
    % Store data in cells %
  	% PATTERN
	PAT.Time        {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = pat.Time;
	PAT.Pos         {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = pat.Pos;
	PAT.Vel         {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = pat.Vel;
   	PAT.VelMean     {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = pat.VelMean;
	PAT.VelSTD      {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = pat.VelSTD;
    PAT.Freq        {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = pat.Freq;
	PAT.Mag         {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = pat.Mag;
	PAT.Phase       {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = pat.Phase;
    % HEAD
	HEAD.Time       {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = head.Time;
	HEAD.Pos        {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = head.Pos;
    HEAD.Vel        {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = head.Vel;
    HEAD.VelMean    {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = head.VelMean;
    HEAD.VelSTD     {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = head.VelSTD;
	HEAD.Err.Pos    {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = head.Err.Pos;
    HEAD.Err.Freq 	{FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = head.Err.Freq;
    HEAD.Err.Mag    {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = head.Err.Mag;
 	HEAD.Err.Phase 	{FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = head.Err.Phase;
	HEAD.Freq       {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = head.Freq;
	HEAD.Mag        {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = head.Mag;
	HEAD.Phase      {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = head.Phase;
	HEAD.COHR.Freq  {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = head.cohr.f;
    HEAD.COHR.Mag   {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = head.cohr.mag;
    HEAD.GAIN       {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = head.GAIN;
	HEAD.PHASE      {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = head.PHASE;
    % WINGS
	WING.Pos        {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = wing.Pos;
	WING.Time       {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = wing.Time;
	WING.Vel        {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = wing.Vel;
	WING.VelMean	{FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = wing.VelMean;
	WING.VelSTD     {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = wing.VelSTD;
	WING.Freq       {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = wing.Freq;
	WING.Mag        {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = wing.Mag;
	WING.Phase      {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = wing.Phase;
	WING.COHR.Freq  {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = wing.cohr.f;
    WING.COHR.Mag   {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = wing.cohr.mag;
	WING.GAIN       {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = wing.GAIN;
  	WING.PHASE      {FD.idxFly(kk),1}{FD.idxAmp(kk),1}(:,end+1) = wing.PHASE;
end
clear jj kk a b t_p t_v hAngles data
disp('DONE')
%% FLY Mean, Median, & STD %%
%---------------------------------------------------------------------------------------------------------------------------------
for kk = 1:n.Fly
    % MEDIAN
    %-----------------------------------------------------------------------------------------------------------------------------
        % PATTERN
        PAT.FlyMed.Time         {kk,1}	= cell2mat(cellfun(@(x) median(x,2),PAT.Time{kk},'UniformOutput',false)');
        PAT.FlyMed.Pos          {kk,1}	= cell2mat(cellfun(@(x) median(x,2),PAT.Pos{kk},'UniformOutput',false)');
        PAT.FlyMed.Vel          {kk,1}	= cell2mat(cellfun(@(x) median(x,2),PAT.Vel{kk},'UniformOutput',false)');
        PAT.FlyMed.Freq         {kk,1}	= cell2mat(cellfun(@(x) median(x,2),PAT.Freq{kk},'UniformOutput',false)');
        PAT.FlyMed.Mag          {kk,1}	= cell2mat(cellfun(@(x) median(x,2),PAT.Mag{kk},'UniformOutput',false)');
        PAT.FlyMed.Phase        {kk,1}	= cell2mat(cellfun(@(x) median(x,2),PAT.Phase{kk},'UniformOutput',false)');
        % HEAD
        HEAD.FlyMed.Time        {kk,1}	= cell2mat(cellfun(@(x) median(x,2),HEAD.Time{kk},'UniformOutput',false)');
        HEAD.FlyMed.Pos         {kk,1}	= cell2mat(cellfun(@(x) median(x,2),HEAD.Pos{kk},'UniformOutput',false)');
        HEAD.FlyMed.Vel         {kk,1}	= cell2mat(cellfun(@(x) median(x,2),HEAD.Vel{kk},'UniformOutput',false)');
        HEAD.FlyMed.Freq        {kk,1}	= cell2mat(cellfun(@(x) median(x,2),HEAD.Freq{kk},'UniformOutput',false)');
        HEAD.FlyMed.Mag         {kk,1}	= cell2mat(cellfun(@(x) median(x,2),HEAD.Mag{kk},'UniformOutput',false)');
        HEAD.FlyMed.Phase       {kk,1}	= cell2mat(cellfun(@(x) median(x,2),HEAD.Phase{kk},'UniformOutput',false)');
        HEAD.FlyMed.COHR.Freq 	{kk,1}	= cell2mat(cellfun(@(x) median(x,2),HEAD.COHR.Freq{kk},'UniformOutput',false)');
        HEAD.FlyMed.COHR.Mag 	{kk,1}	= cell2mat(cellfun(@(x) median(x,2),HEAD.COHR.Mag{kk},'UniformOutput',false)');
        HEAD.FlyMed.GAIN        {kk,1}	= cell2mat(cellfun(@(x) median(x,2),HEAD.Phase{kk},'UniformOutput',false)');
        HEAD.FlyMed.PHASE    	{kk,1}	= cell2mat(cellfun(@(x) median(x,2),HEAD.PHASE{kk},'UniformOutput',false)');

        HEAD.FlyMed.Err.Pos    	{kk,1}	= cell2mat(cellfun(@(x) median(x,2),HEAD.Err.Pos{kk},'UniformOutput',false)');
        HEAD.FlyMed.Err.Freq  	{kk,1}	= cell2mat(cellfun(@(x) median(x,2),HEAD.Err.Freq{kk},'UniformOutput',false)');
        HEAD.FlyMed.Err.Mag    	{kk,1}	= cell2mat(cellfun(@(x) median(x,2),HEAD.Err.Mag{kk},'UniformOutput',false)');
        HEAD.FlyMed.Err.Phase  	{kk,1}	= cell2mat(cellfun(@(x) median(x,2),HEAD.Err.Phase{kk},'UniformOutput',false)');
        % WINGS
        WING.FlyMed.Time        {kk,1}	= cell2mat(cellfun(@(x) median(x,2),WING.Time{kk},'UniformOutput',false)');
        WING.FlyMed.Pos         {kk,1}	= cell2mat(cellfun(@(x) median(x,2),WING.Pos{kk},'UniformOutput',false)');
        WING.FlyMed.Vel         {kk,1}	= cell2mat(cellfun(@(x) median(x,2),WING.Vel{kk},'UniformOutput',false)');
        WING.FlyMed.Freq        {kk,1}	= cell2mat(cellfun(@(x) median(x,2),WING.Freq{kk},'UniformOutput',false)');
        WING.FlyMed.Mag         {kk,1}	= cell2mat(cellfun(@(x) median(x,2),WING.Mag{kk},'UniformOutput',false)');
        WING.FlyMed.Phase       {kk,1}	= cell2mat(cellfun(@(x) median(x,2),WING.Phase{kk},'UniformOutput',false)');
        WING.FlyMed.COHR.Freq 	{kk,1}	= cell2mat(cellfun(@(x) median(x,2),WING.COHR.Freq{kk},'UniformOutput',false)');
        WING.FlyMed.COHR.Mag 	{kk,1}	= cell2mat(cellfun(@(x) median(x,2),WING.COHR.Mag{kk},'UniformOutput',false)');
        WING.FlyMed.GAIN        {kk,1}	= cell2mat(cellfun(@(x) median(x,2),WING.Phase{kk},'UniformOutput',false)');
        WING.FlyMed.PHASE    	{kk,1}	= cell2mat(cellfun(@(x) median(x,2),WING.PHASE{kk},'UniformOutput',false)');
    
    % MEAN
    %-----------------------------------------------------------------------------------------------------------------------------
        % PATTERN
        PAT.FlyMean.Time    	{kk,1}	= cell2mat(cellfun(@(x) mean(x,2),PAT.Time{kk},'UniformOutput',false)');
        PAT.FlyMean.Pos       	{kk,1}	= cell2mat(cellfun(@(x) mean(x,2),PAT.Pos{kk},'UniformOutput',false)');
        PAT.FlyMean.Vel       	{kk,1}	= cell2mat(cellfun(@(x) mean(x,2),PAT.Vel{kk},'UniformOutput',false)');
        PAT.FlyMean.Freq     	{kk,1}	= cell2mat(cellfun(@(x) mean(x,2),PAT.Freq{kk},'UniformOutput',false)');
        PAT.FlyMean.Mag      	{kk,1}	= cell2mat(cellfun(@(x) mean(x,2),PAT.Mag{kk},'UniformOutput',false)');
        PAT.FlyMean.Phase     	{kk,1}	= cell2mat(cellfun(@(x) mean(x,2),PAT.Phase{kk},'UniformOutput',false)');
        % HEAD
        HEAD.FlyMean.Time     	{kk,1}	= cell2mat(cellfun(@(x) mean(x,2),HEAD.Time{kk},'UniformOutput',false)');
        HEAD.FlyMean.Pos      	{kk,1}	= cell2mat(cellfun(@(x) mean(x,2),HEAD.Pos{kk},'UniformOutput',false)');
        HEAD.FlyMean.Vel     	{kk,1}	= cell2mat(cellfun(@(x) mean(x,2),HEAD.Vel{kk},'UniformOutput',false)');
        HEAD.FlyMean.Freq    	{kk,1}	= cell2mat(cellfun(@(x) mean(x,2),HEAD.Freq{kk},'UniformOutput',false)');
        HEAD.FlyMean.Mag      	{kk,1}	= cell2mat(cellfun(@(x) mean(x,2),HEAD.Mag{kk},'UniformOutput',false)');
        HEAD.FlyMean.Phase  	{kk,1}	= cell2mat(cellfun(@(x) mean(x,2),HEAD.Phase{kk},'UniformOutput',false)');
        HEAD.FlyMean.COHR.Freq 	{kk,1}	= cell2mat(cellfun(@(x) mean(x,2),HEAD.COHR.Freq{kk},'UniformOutput',false)');
        HEAD.FlyMean.COHR.Mag 	{kk,1}	= cell2mat(cellfun(@(x) mean(x,2),HEAD.COHR.Mag{kk},'UniformOutput',false)');
        HEAD.FlyMean.GAIN    	{kk,1}	= cell2mat(cellfun(@(x) mean(x,2),HEAD.Phase{kk},'UniformOutput',false)');
        HEAD.FlyMean.PHASE   	{kk,1}	= cell2mat(cellfun(@(x) mean(x,2),HEAD.PHASE{kk},'UniformOutput',false)');

        HEAD.FlyMean.Err.Pos  	{kk,1}	= cell2mat(cellfun(@(x) mean(x,2),HEAD.Err.Pos{kk},'UniformOutput',false)');
        HEAD.FlyMean.Err.Freq  	{kk,1}	= cell2mat(cellfun(@(x) mean(x,2),HEAD.Err.Freq{kk},'UniformOutput',false)');
        HEAD.FlyMean.Err.Mag 	{kk,1}	= cell2mat(cellfun(@(x) mean(x,2),HEAD.Err.Mag{kk},'UniformOutput',false)');
        HEAD.FlyMean.Err.Phase	{kk,1}	= cell2mat(cellfun(@(x) mean(x,2),HEAD.Err.Phase{kk},'UniformOutput',false)');
        % WINGS
        WING.FlyMean.Time      	{kk,1}	= cell2mat(cellfun(@(x) mean(x,2),WING.Time{kk},'UniformOutput',false)');
        WING.FlyMean.Pos    	{kk,1}	= cell2mat(cellfun(@(x) mean(x,2),WING.Pos{kk},'UniformOutput',false)');
        WING.FlyMean.Vel      	{kk,1}	= cell2mat(cellfun(@(x) mean(x,2),WING.Vel{kk},'UniformOutput',false)');
        WING.FlyMean.Freq   	{kk,1}	= cell2mat(cellfun(@(x) mean(x,2),WING.Freq{kk},'UniformOutput',false)');
        WING.FlyMean.Mag     	{kk,1}	= cell2mat(cellfun(@(x) mean(x,2),WING.Mag{kk},'UniformOutput',false)');
        WING.FlyMean.Phase   	{kk,1}	= cell2mat(cellfun(@(x) mean(x,2),WING.Phase{kk},'UniformOutput',false)');
        WING.FlyMean.COHR.Freq	{kk,1}	= cell2mat(cellfun(@(x) mean(x,2),WING.COHR.Freq{kk},'UniformOutput',false)');
        
	% STD
    %-----------------------------------------------------------------------------------------------------------------------------
        % PATTERN
        PAT.FlySTD.Time         {kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),PAT.Time{kk},'UniformOutput',false)');
        PAT.FlySTD.Pos       	{kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),PAT.Pos{kk},'UniformOutput',false)');
        PAT.FlySTD.Vel       	{kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),PAT.Vel{kk},'UniformOutput',false)');
        PAT.FlySTD.Freq     	{kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),PAT.Freq{kk},'UniformOutput',false)');
        PAT.FlySTD.Mag      	{kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),PAT.Mag{kk},'UniformOutput',false)');
        PAT.FlySTD.Phase     	{kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),PAT.Phase{kk},'UniformOutput',false)');
        % HEAD
        HEAD.FlySTD.Time     	{kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),HEAD.Time{kk},'UniformOutput',false)');
        HEAD.FlySTD.Pos      	{kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),HEAD.Pos{kk},'UniformOutput',false)');
        HEAD.FlySTD.Vel     	{kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),HEAD.Vel{kk},'UniformOutput',false)');
        HEAD.FlySTD.Freq    	{kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),HEAD.Freq{kk},'UniformOutput',false)');
        HEAD.FlySTD.Mag      	{kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),HEAD.Mag{kk},'UniformOutput',false)');
        HEAD.FlySTD.Phase       {kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),HEAD.Phase{kk},'UniformOutput',false)');
        HEAD.FlySTD.COHR.Freq 	{kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),HEAD.COHR.Freq{kk},'UniformOutput',false)');
        HEAD.FlySTD.COHR.Mag 	{kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),HEAD.COHR.Mag{kk},'UniformOutput',false)');
        HEAD.FlySTD.GAIN    	{kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),HEAD.Phase{kk},'UniformOutput',false)');
        HEAD.FlySTD.PHASE   	{kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),HEAD.PHASE{kk},'UniformOutput',false)');

        HEAD.FlySTD.Err.Pos  	{kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),HEAD.Err.Pos{kk},'UniformOutput',false)');
        HEAD.FlySTD.Err.Freq  	{kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),HEAD.Err.Freq{kk},'UniformOutput',false)');
        HEAD.FlySTD.Err.Mag 	{kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),HEAD.Err.Mag{kk},'UniformOutput',false)');
        HEAD.FlySTD.Err.Phase	{kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),HEAD.Err.Phase{kk},'UniformOutput',false)');
        % WINGS
        WING.FlySTD.Time      	{kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),WING.Time{kk},'UniformOutput',false)');
        WING.FlySTD.Pos         {kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),WING.Pos{kk},'UniformOutput',false)');
        WING.FlySTD.Vel      	{kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),WING.Vel{kk},'UniformOutput',false)');
        WING.FlySTD.Freq        {kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),WING.Freq{kk},'UniformOutput',false)');
        WING.FlySTD.Mag     	{kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),WING.Mag{kk},'UniformOutput',false)');
        WING.FlySTD.Phase   	{kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),WING.Phase{kk},'UniformOutput',false)');
        WING.FlySTD.COHR.Freq	{kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),WING.COHR.Freq{kk},'UniformOutput',false)');
        WING.FlySTD.COHR.Mag 	{kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),WING.COHR.Mag{kk},'UniformOutput',false)');
        WING.FlySTD.GAIN     	{kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),WING.Phase{kk},'UniformOutput',false)');
        WING.FlySTD.PHASE    	{kk,1}	= cell2mat(cellfun(@(x) std(x,0,3),WING.PHASE{kk},'UniformOutput',false)');
end
clear kk
%% GRAND Mean, Median, & STD %%
%---------------------------------------------------------------------------------------------------------------------------------
% MEDIAN
%---------------------------------------------------------------------------------------------------------------------------------
    PAT.GrandMed.Time           = median(cat(3,PAT.FlyMed.Time{:}),3);
    PAT.GrandMed.Pos            = median(cat(3,PAT.FlyMed.Pos{:}),3);
    PAT.GrandMed.Vel            = median(cat(3,PAT.FlyMed.Vel{:}),3);
    PAT.GrandMed.Freq           = median(cat(3,PAT.FlyMed.Freq{:}),3);
    PAT.GrandMed.Phase          = median(cat(3,PAT.FlyMed.Phase{:}),3);
    PAT.GrandMed.Mag            = median(cat(3,PAT.FlyMed.Mag{:}),3);

    HEAD.GrandMed.Time          = median(cat(3,HEAD.FlyMed.Time{:}),3);
    HEAD.GrandMed.Pos           = median(cat(3,HEAD.FlyMed.Pos{:}),3);
    HEAD.GrandMed.Vel           = median(cat(3,HEAD.FlyMed.Vel{:}),3);
    HEAD.GrandMed.Freq          = median(cat(3,HEAD.FlyMed.Freq{:}),3);
    HEAD.GrandMed.Phase         = median(cat(3,HEAD.FlyMed.Phase{:}),3);
    HEAD.GrandMed.Mag           = median(cat(3,HEAD.FlyMed.Mag{:}),3);
    HEAD.GrandMed.COHR.Freq   	= median(cat(3,HEAD.FlyMed.COHR.Freq{:}),3);
    HEAD.GrandMed.COHR.Mag   	= median(cat(3,HEAD.FlyMed.COHR.Mag{:}),3);
    HEAD.GrandMed.GAIN          = median(cat(3,HEAD.FlyMed.GAIN{:}),3);
    HEAD.GrandMed.PHASE         = median(cat(3,HEAD.FlyMed.PHASE{:}),3);
    HEAD.GrandMed.Err.Pos     	= median(cat(3,HEAD.FlyMed.Err.Pos {:}),3);
    HEAD.GrandMed.Err.Freq     	= median(cat(3,HEAD.FlyMed.Err.Freq {:}),3);
    HEAD.GrandMed.Err.Mag     	= median(cat(3,HEAD.FlyMed.Err.Mag {:}),3);
    HEAD.GrandMed.Err.Phase    	= median(cat(3,HEAD.FlyMed.Err.Phase {:}),3);

    WING.GrandMed.Time          = median(cat(3,WING.FlyMed.Time{:}),3);
    WING.GrandMed.Pos           = median(cat(3,WING.FlyMed.Pos{:}),3);
    WING.GrandMed.Vel           = median(cat(3,WING.FlyMed.Vel{:}),3);
    WING.GrandMed.Freq          = median(cat(3,WING.FlyMed.Freq{:}),3);
    WING.GrandMed.Phase         = median(cat(3,WING.FlyMed.Phase{:}),3);
    WING.GrandMed.Mag           = median(cat(3,WING.FlyMed.Mag{:}),3);
    WING.GrandMed.COHR.Freq   	= median(cat(3,WING.FlyMed.COHR.Freq{:}),3);
    WING.GrandMed.COHR.Mag   	= median(cat(3,WING.FlyMed.COHR.Mag{:}),3);
    WING.GrandMed.GAIN          = median(cat(3,WING.FlyMed.GAIN{:}),3);
    WING.GrandMed.PHASE         = median(cat(3,WING.FlyMed.PHASE{:}),3);

% MEAN
%---------------------------------------------------------------------------------------------------------------------------------
    PAT.GrandMed.Time           = mean(cat(3,PAT.FlyMed.Time{:}),3);
    PAT.GrandMed.Pos            = mean(cat(3,PAT.FlyMed.Pos{:}),3);
    PAT.GrandMed.Vel            = mean(cat(3,PAT.FlyMed.Vel{:}),3);
    PAT.GrandMed.Freq           = mean(cat(3,PAT.FlyMed.Freq{:}),3);
    PAT.GrandMed.Phase          = mean(cat(3,PAT.FlyMed.Phase{:}),3);
    PAT.GrandMed.Mag            = mean(cat(3,PAT.FlyMed.Mag{:}),3);

    HEAD.GrandMed.Time          = mean(cat(3,HEAD.FlyMed.Time{:}),3);
    HEAD.GrandMed.Pos           = mean(cat(3,HEAD.FlyMed.Pos{:}),3);
    HEAD.GrandMed.Vel           = mean(cat(3,HEAD.FlyMed.Vel{:}),3);
    HEAD.GrandMed.Freq          = mean(cat(3,HEAD.FlyMed.Freq{:}),3);
    HEAD.GrandMed.Phase         = mean(cat(3,HEAD.FlyMed.Phase{:}),3);
    HEAD.GrandMed.Mag           = mean(cat(3,HEAD.FlyMed.Mag{:}),3);
    HEAD.GrandMed.COHR.Freq   	= mean(cat(3,HEAD.FlyMed.COHR.Freq{:}),3);
    HEAD.GrandMed.COHR.Mag   	= mean(cat(3,HEAD.FlyMed.COHR.Mag{:}),3);
    HEAD.GrandMed.GAIN          = mean(cat(3,HEAD.FlyMed.GAIN{:}),3);
    HEAD.GrandMed.PHASE         = mean(cat(3,HEAD.FlyMed.PHASE{:}),3);
    HEAD.GrandMed.Err.Pos     	= mean(cat(3,HEAD.FlyMed.Err.Pos {:}),3);
    HEAD.GrandMed.Err.Freq     	= mean(cat(3,HEAD.FlyMed.Err.Freq {:}),3);
    HEAD.GrandMed.Err.Mag     	= mean(cat(3,HEAD.FlyMed.Err.Mag {:}),3);
    HEAD.GrandMed.Err.Phase    	= mean(cat(3,HEAD.FlyMed.Err.Phase {:}),3);

    WING.GrandMed.Time          = mean(cat(3,WING.FlyMed.Time{:}),3);
    WING.GrandMed.Pos           = mean(cat(3,WING.FlyMed.Pos{:}),3);
    WING.GrandMed.Vel           = mean(cat(3,WING.FlyMed.Vel{:}),3);
    WING.GrandMed.Freq          = mean(cat(3,WING.FlyMed.Freq{:}),3);
    WING.GrandMed.Phase         = mean(cat(3,WING.FlyMed.Phase{:}),3);
    WING.GrandMed.Mag           = mean(cat(3,WING.FlyMed.Mag{:}),3);
    WING.GrandMed.COHR.Freq   	= mean(cat(3,WING.FlyMed.COHR.Freq{:}),3);
    WING.GrandMed.COHR.Mag   	= mean(cat(3,WING.FlyMed.COHR.Mag{:}),3);
    WING.GrandMed.GAIN          = mean(cat(3,WING.FlyMed.GAIN{:}),3);
    WING.GrandMed.PHASE         = mean(cat(3,WING.FlyMed.PHASE{:}),3);

% STD
%---------------------------------------------------------------------------------------------------------------------------------
    PAT.GrandSTD.Time           = std(cat(3,PAT.FlyMed.Time{:}),0,3);
    PAT.GrandSTD.Pos            = std(cat(3,PAT.FlyMed.Pos{:}),0,3);
    PAT.GrandSTD.Vel            = std(cat(3,PAT.FlyMed.Vel{:}),0,3);
    PAT.GrandSTD.Freq           = std(cat(3,PAT.FlyMed.Freq{:}),0,3);
    PAT.GrandSTD.Phase          = std(cat(3,PAT.FlyMed.Phase{:}),0,3);
    PAT.GrandSTD.Mag            = std(cat(3,PAT.FlyMed.Mag{:}),0,3);

    HEAD.GrandSTD.Time          = std(cat(3,HEAD.FlyMed.Time{:}),0,3);
    HEAD.GrandSTD.Pos           = std(cat(3,HEAD.FlyMed.Pos{:}),0,3);
    HEAD.GrandSTD.Vel           = std(cat(3,HEAD.FlyMed.Vel{:}),0,3);
    HEAD.GrandSTD.Freq          = std(cat(3,HEAD.FlyMed.Freq{:}),0,3);
    HEAD.GrandSTD.Phase         = std(cat(3,HEAD.FlyMed.Phase{:}),0,3);
    HEAD.GrandSTD.Mag           = std(cat(3,HEAD.FlyMed.Mag{:}),0,3);
    HEAD.GrandSTD.COHR.Freq   	= std(cat(3,HEAD.FlyMed.COHR.Freq{:}),0,3);
    HEAD.GrandSTD.COHR.Mag   	= std(cat(3,HEAD.FlyMed.COHR.Mag{:}),0,3);
    HEAD.GrandSTD.GAIN          = std(cat(3,HEAD.FlyMed.GAIN{:}),0,3);
    HEAD.GrandSTD.PHASE         = std(cat(3,HEAD.FlyMed.PHASE{:}),0,3);
    HEAD.GrandSTD.Err.Pos     	= std(cat(3,HEAD.FlyMed.Err.Pos {:}),0,3);
    HEAD.GrandSTD.Err.Freq     	= std(cat(3,HEAD.FlyMed.Err.Freq {:}),0,3);
    HEAD.GrandSTD.Err.Mag     	= std(cat(3,HEAD.FlyMed.Err.Mag {:}),0,3);
    HEAD.GrandSTD.Err.Phase    	= std(cat(3,HEAD.FlyMed.Err.Phase {:}),0,3);

    WING.GrandSTD.Time          = std(cat(3,WING.FlyMed.Time{:}),0,3);
    WING.GrandSTD.Pos           = std(cat(3,WING.FlyMed.Pos{:}),0,3);
    WING.GrandSTD.Vel           = std(cat(3,WING.FlyMed.Vel{:}),0,3);
    WING.GrandSTD.Freq          = std(cat(3,WING.FlyMed.Freq{:}),0,3);
    WING.GrandSTD.Phase         = std(cat(3,WING.FlyMed.Phase{:}),0,3);
    WING.GrandSTD.Mag           = std(cat(3,WING.FlyMed.Mag{:}),0,3);
    WING.GrandSTD.COHR.Freq   	= std(cat(3,WING.FlyMed.COHR.Freq{:}),0,3);
    WING.GrandSTD.COHR.Mag   	= std(cat(3,WING.FlyMed.COHR.Mag{:}),0,3);
    WING.GrandSTD.GAIN          = std(cat(3,WING.FlyMed.GAIN{:}),0,3);
    WING.GrandSTD.PHASE         = std(cat(3,WING.FlyMed.PHASE{:}),0,3);
%% Save ouputs as structure %%
%---------------------------------------------------------------------------------------------------------------------------------
save([rootdir 'DATA\DATA.mat'],'PAT','WING','HEAD','FD','T','n','unq')
end
