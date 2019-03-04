function [PAT,WING,FD,T,n,unq] = MakeData_Chirp_HeadFixed(rootdir,filename)
%% MakeData_Chirp_HeadFree: Reads in all raw trials, transforms data, and saves in organized structure for use with figure functions
%   INPUTS:
%       root    : root directory
%   OUTPUTS:
%       PAT     : pattern structure
%       WING   	: wings structure
%       FD      : file data
%       T       : fly data table
%       n       : field #'s 
%       unq   	: unique fields
%---------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE INPUT %
% rootdir = 'H:\Experiment_HeadExcitation\Chirp\HeadFixed\';
% filename = 'Chirp_HeadFixed_DATA';
%---------------------------------------------------------------------------------------------------------------------------------
%% Setup Directories %%
%---------------------------------------------------------------------------------------------------------------------------------
root.daq = rootdir;

% Select files
[FILES, PATH.daq] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select DAQ trials', root.daq, 'MultiSelect','on');
FILES = FILES';

%% Process File Data %%
%---------------------------------------------------------------------------------------------------------------------------------
% Read in data from DAQ file names
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
disp('Loading...')
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
    end
end
% Store data in organized cells
tt = (0:1/200:20)';
tt = tt(1:end-1);
for kk = 1:n.Trial
    disp(kk)
    % Load DAQ data %
    data = [];
	load([PATH.daq   FILES{kk}],'data','t_p'); % load pattern x-position
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
 	pat.Pos         = FitPanel(pat.Pos,pat.Time,tt); % fit panel data
    pat.Time        = tt; % set new panel time
	pat.Vel         = [diff(pat.Pos)/(1/pat.Fs) ; 0]; % pattern velocity [deg/s]
    pat.VelMean     = mean(abs(pat.Vel)); % mean pattern velocity [deg/s]
	pat.VelSTD      = mean(abs(pat.Vel)); % STD pattern velocity [deg/s]
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
	wing.Time       = tt;
	wing.Pos        = resample(wing.Pos,4000,wing.n);
	wing.Vel        = resample(wing.Vel,4000,wing.n);
    wing.Fs         = 200;
    %-----------------------------------------------------------------------------------------------------------------------------
    % Convert wings & pattern data to frequency domain using FFT %
    [wing.Freq , wing.Mag , wing.Phase]     = FFT(wing.Time,wing.Pos);
    [pat.Freq  , pat.Mag  , pat.Phase ] 	= FFT(pat.Time,pat.Pos);
    %-----------------------------------------------------------------------------------------------------------------------------
    % Calculate coherence %
	[wing.cohr.mag,wing.cohr.f] = mscohere(pat.Pos , wing.Pos ,[],[] , wing.Freq , wing.Fs);
    %-----------------------------------------------------------------------------------------------------------------------------
	% Calculate BODE gain & phase difference for wings %
 	wing.GAIN   = medfilt1(wing.Mag./pat.Mag,10);
    wing.PHASE  = medfilt1(-(pat.Phase - wing.Phase),10);
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
%% FLY Stats by Fly %%
%---------------------------------------------------------------------------------------------------------------------------------
for kk = 1:n.Fly
    for jj = 1:n.Amp
    % MEDIAN
    %-----------------------------------------------------------------------------------------------------------------------------
        % PATTERN
        PAT.FlyMed.Time         {kk,1}(:,jj)	= median(PAT.Time{kk}{jj},2);
        PAT.FlyMed.Pos          {kk,1}(:,jj)	= median(PAT.Pos{kk}{jj},2);
        PAT.FlyMed.Vel          {kk,1}(:,jj)	= median(PAT.Vel{kk}{jj},2);
        PAT.FlyMed.Freq         {kk,1}(:,jj)	= median(PAT.Freq{kk}{jj},2);
        PAT.FlyMed.Mag          {kk,1}(:,jj)	= median(PAT.Mag{kk}{jj},2);
        PAT.FlyMed.Phase     	{kk,1}(:,jj)	= circ_median(PAT.Phase{kk}{jj},2)';
        % WINGS
        WING.FlyMed.Time    	{kk,1}(:,jj)	= median(WING.Time{kk}{jj},2);
        WING.FlyMed.Pos         {kk,1}(:,jj)	= median(WING.Pos{kk}{jj},2);
        WING.FlyMed.Vel         {kk,1}(:,jj)	= median(WING.Vel{kk}{jj},2);
        WING.FlyMed.Freq    	{kk,1}(:,jj)	= median(WING.Freq{kk}{jj},2);
        WING.FlyMed.Mag     	{kk,1}(:,jj)	= median(WING.Mag{kk}{jj},2);
        WING.FlyMed.Phase    	{kk,1}(:,jj)	= circ_median(WING.Phase{kk}{jj},2)';
        WING.FlyMed.COHR.Freq 	{kk,1}(:,jj)	= median(WING.COHR.Freq{kk}{jj},2);
        WING.FlyMed.COHR.Mag 	{kk,1}(:,jj)	= median(WING.COHR.Mag{kk}{jj},2);
        WING.FlyMed.GAIN     	{kk,1}(:,jj)	= median(WING.GAIN{kk}{jj},2);
        WING.FlyMed.PHASE    	{kk,1}(:,jj)	= circ_median(WING.PHASE{kk}{jj},2)';
        
	% MEAN
    %-----------------------------------------------------------------------------------------------------------------------------
        % PATTERN
        PAT.FlyMean.Time    	{kk,1}(:,jj)	= mean(PAT.Time{kk}{jj},2);
        PAT.FlyMean.Pos       	{kk,1}(:,jj)	= mean(PAT.Pos{kk}{jj},2);
        PAT.FlyMean.Vel      	{kk,1}(:,jj)	= mean(PAT.Vel{kk}{jj},2);
        PAT.FlyMean.Freq     	{kk,1}(:,jj)	= mean(PAT.Freq{kk}{jj},2);
        PAT.FlyMean.Mag       	{kk,1}(:,jj)	= mean(PAT.Mag{kk}{jj},2);
        PAT.FlyMean.Phase     	{kk,1}(:,jj)	= circ_mean(PAT.Phase{kk}{jj},[],2)';
        % WINGS
        WING.FlyMean.Time    	{kk,1}(:,jj)	= mean(WING.Time{kk}{jj},2);
        WING.FlyMean.Pos      	{kk,1}(:,jj)	= mean(WING.Pos{kk}{jj},2);
        WING.FlyMean.Vel       	{kk,1}(:,jj)	= mean(WING.Vel{kk}{jj},2);
        WING.FlyMean.Freq    	{kk,1}(:,jj)	= mean(WING.Freq{kk}{jj},2);
        WING.FlyMean.Mag     	{kk,1}(:,jj)	= mean(WING.Mag{kk}{jj},2);
        WING.FlyMean.Phase    	{kk,1}(:,jj)	= circ_mean(WING.Phase{kk}{jj},[],2)';
        WING.FlyMean.COHR.Freq 	{kk,1}(:,jj)	= mean(WING.COHR.Freq{kk}{jj},2);
        WING.FlyMean.COHR.Mag 	{kk,1}(:,jj)	= mean(WING.COHR.Mag{kk}{jj},2);
        WING.FlyMean.GAIN     	{kk,1}(:,jj)	= mean(WING.GAIN{kk}{jj},2);
        WING.FlyMean.PHASE    	{kk,1}(:,jj)	= circ_mean(WING.PHASE{kk}{jj},[],2)';
    % STD
    %-----------------------------------------------------------------------------------------------------------------------------
        % PATTERN
        PAT.FlySTD.Time         {kk,1}(:,jj)	= std(PAT.Time{kk}{jj},0,2);
        PAT.FlySTD.Pos          {kk,1}(:,jj)	= std(PAT.Pos{kk}{jj},0,2);
        PAT.FlySTD.Vel          {kk,1}(:,jj)	= std(PAT.Vel{kk}{jj},0,2);
        PAT.FlySTD.Freq         {kk,1}(:,jj)	= std(PAT.Freq{kk}{jj},0,2);
        PAT.FlySTD.Mag          {kk,1}(:,jj)	= std(PAT.Mag{kk}{jj},0,2);
        PAT.FlySTD.Phase     	{kk,1}(:,jj)	= circ_std(PAT.Phase{kk}{jj},[],[],2);
        % WINGS
        WING.FlySTD.Time    	{kk,1}(:,jj)	= std(WING.Time{kk}{jj},0,2);
        WING.FlySTD.Pos         {kk,1}(:,jj)	= std(WING.Pos{kk}{jj},0,2);
        WING.FlySTD.Vel         {kk,1}(:,jj)	= std(WING.Vel{kk}{jj},0,2);
        WING.FlySTD.Freq    	{kk,1}(:,jj)	= std(WING.Freq{kk}{jj},0,2);
        WING.FlySTD.Mag     	{kk,1}(:,jj)	= std(WING.Mag{kk}{jj},0,2);
        WING.FlySTD.Phase    	{kk,1}(:,jj)	= circ_std(WING.Phase{kk}{jj},[],[],2);
        WING.FlySTD.COHR.Freq 	{kk,1}(:,jj)	= std(WING.COHR.Freq{kk}{jj},0,2);
        WING.FlySTD.COHR.Mag 	{kk,1}(:,jj)	= std(WING.COHR.Mag{kk}{jj},0,2);
        WING.FlySTD.GAIN     	{kk,1}(:,jj)	= std(WING.GAIN{kk}{jj},0,2);
        WING.FlySTD.PHASE    	{kk,1}(:,jj)	= circ_std(WING.PHASE{kk}{jj},[],[],2);
    end
end
clear kk
%% FLY stats by Amplitude %%
%---------------------------------------------------------------------------------------------------------------------------------
for jj = 1:n.Amp
    for kk = 1:n.Fly
	% MEDIAN
    %-----------------------------------------------------------------------------------------------------------------------------
        PAT.FlyAmpMed.Time          {jj,1}(:,kk)  	= PAT.FlyMed.Time{kk}(:,jj);
        PAT.FlyAmpMed.Pos           {jj,1}(:,kk) 	= PAT.FlyMed.Pos{kk}(:,jj);
    	PAT.FlyAmpMed.Vel           {jj,1}(:,kk)  	= PAT.FlyMed.Vel{kk}(:,jj);
        PAT.FlyAmpMed.Freq          {jj,1}(:,kk)  	= PAT.FlyMed.Freq{kk}(:,jj);
        PAT.FlyAmpMed.Phase         {jj,1}(:,kk)  	= PAT.FlyMed.Phase{kk}(:,jj);
        PAT.FlyAmpMed.Mag           {jj,1}(:,kk)   	= PAT.FlyMed.Mag{kk}(:,jj);

        WING.FlyAmpMed.Time         {jj,1}(:,kk)  	= WING.FlyMed.Time{kk}(:,jj);
        WING.FlyAmpMed.Pos          {jj,1}(:,kk)   	= WING.FlyMed.Pos{kk}(:,jj);
    	WING.FlyAmpMed.Vel          {jj,1}(:,kk)   	= WING.FlyMed.Vel{kk}(:,jj);
        WING.FlyAmpMed.Freq         {jj,1}(:,kk)   	= WING.FlyMed.Freq{kk}(:,jj);
        WING.FlyAmpMed.Phase        {jj,1}(:,kk)  	= WING.FlyMed.Phase{kk}(:,jj);
        WING.FlyAmpMed.Mag          {jj,1}(:,kk)  	= WING.FlyMed.Mag{kk}(:,jj);
        WING.FlyAmpMed.COHR.Freq    {jj,1}(:,kk)	= WING.FlyMed.COHR.Freq{kk}(:,jj);
        WING.FlyAmpMed.COHR.Mag     {jj,1}(:,kk)	= WING.FlyMed.COHR.Mag{kk}(:,jj);
        WING.FlyAmpMed.GAIN         {jj,1}(:,kk)   	= WING.FlyMed.GAIN{kk}(:,jj);
        WING.FlyAmpMed.PHASE        {jj,1}(:,kk)   	= WING.FlyMed.PHASE{kk}(:,jj);
	% MEAN
    %-----------------------------------------------------------------------------------------------------------------------------
        PAT.FlyAmpMean.Time        	{jj,1}(:,kk)  	= PAT.FlyMean.Time{kk}(:,jj);
        PAT.FlyAmpMean.Pos         	{jj,1}(:,kk) 	= PAT.FlyMean.Pos{kk}(:,jj);
    	PAT.FlyAmpMean.Vel         	{jj,1}(:,kk)  	= PAT.FlyMean.Vel{kk}(:,jj);
        PAT.FlyAmpMean.Freq      	{jj,1}(:,kk)  	= PAT.FlyMean.Freq{kk}(:,jj);
        PAT.FlyAmpMean.Phase     	{jj,1}(:,kk)  	= PAT.FlyMean.Phase{kk}(:,jj);
        PAT.FlyAmpMean.Mag       	{jj,1}(:,kk)   	= PAT.FlyMean.Mag{kk}(:,jj);
        
        WING.FlyAmpMean.Time     	{jj,1}(:,kk)  	= WING.FlyMean.Time{kk}(:,jj);
        WING.FlyAmpMean.Pos      	{jj,1}(:,kk)   	= WING.FlyMean.Pos{kk}(:,jj);
    	WING.FlyAmpMean.Vel        	{jj,1}(:,kk)   	= WING.FlyMean.Vel{kk}(:,jj);
        WING.FlyAmpMean.Freq    	{jj,1}(:,kk)   	= WING.FlyMean.Freq{kk}(:,jj);
        WING.FlyAmpMean.Phase    	{jj,1}(:,kk)  	= WING.FlyMean.Phase{kk}(:,jj);
        WING.FlyAmpMean.Mag       	{jj,1}(:,kk)  	= WING.FlyMean.Mag{kk}(:,jj);
        WING.FlyAmpMean.COHR.Freq  	{jj,1}(:,kk)	= WING.FlyMean.COHR.Freq{kk}(:,jj);
        WING.FlyAmpMean.COHR.Mag  	{jj,1}(:,kk)	= WING.FlyMean.COHR.Mag{kk}(:,jj);
        WING.FlyAmpMean.GAIN      	{jj,1}(:,kk)   	= WING.FlyMean.GAIN{kk}(:,jj);
        WING.FlyAmpMean.PHASE    	{jj,1}(:,kk)   	= WING.FlyMean.PHASE{kk}(:,jj);
	% STD
    %-----------------------------------------------------------------------------------------------------------------------------
        PAT.FlyAmpSTD.Time        	{jj,1}(:,kk)  	= PAT.FlySTD.Time{kk}(:,jj);
        PAT.FlyAmpSTD.Pos         	{jj,1}(:,kk) 	= PAT.FlySTD.Pos{kk}(:,jj);
    	PAT.FlyAmpSTD.Vel         	{jj,1}(:,kk)  	= PAT.FlySTD.Vel{kk}(:,jj);
        PAT.FlyAmpSTD.Freq      	{jj,1}(:,kk)  	= PAT.FlySTD.Freq{kk}(:,jj);
        PAT.FlyAmpSTD.Phase     	{jj,1}(:,kk)  	= PAT.FlySTD.Phase{kk}(:,jj);
        PAT.FlyAmpSTD.Mag       	{jj,1}(:,kk)   	= PAT.FlySTD.Mag{kk}(:,jj);

        WING.FlyAmpSTD.Time     	{jj,1}(:,kk)  	= WING.FlySTD.Time{kk}(:,jj);
        WING.FlyAmpSTD.Pos      	{jj,1}(:,kk)   	= WING.FlySTD.Pos{kk}(:,jj);
    	WING.FlyAmpSTD.Vel        	{jj,1}(:,kk)   	= WING.FlySTD.Vel{kk}(:,jj);
        WING.FlyAmpSTD.Freq         {jj,1}(:,kk)   	= WING.FlySTD.Freq{kk}(:,jj);
        WING.FlyAmpSTD.Phase    	{jj,1}(:,kk)  	= WING.FlySTD.Phase{kk}(:,jj);
        WING.FlyAmpSTD.Mag       	{jj,1}(:,kk)  	= WING.FlySTD.Mag{kk}(:,jj);
        WING.FlyAmpSTD.COHR.Freq  	{jj,1}(:,kk)	= WING.FlySTD.COHR.Freq{kk}(:,jj);
        WING.FlyAmpSTD.COHR.Mag  	{jj,1}(:,kk)	= WING.FlySTD.COHR.Mag{kk}(:,jj);
        WING.FlyAmpSTD.GAIN      	{jj,1}(:,kk)   	= WING.FlySTD.GAIN{kk}(:,jj);
        WING.FlyAmpSTD.PHASE    	{jj,1}(:,kk)   	= WING.FlySTD.PHASE{kk}(:,jj);
    end
end
%% GRAND Stats %%
%---------------------------------------------------------------------------------------------------------------------------------
% MEDIAN
%---------------------------------------------------------------------------------------------------------------------------------
    PAT.GrandMed.Time           = cell2mat((cellfun(@(x) median(x,2),PAT.FlyAmpMed.Time,'UniformOutput',false))');
    PAT.GrandMed.Pos            = cell2mat((cellfun(@(x) median(x,2),PAT.FlyAmpMed.Pos,'UniformOutput',false))');
    PAT.GrandMed.Vel            = cell2mat((cellfun(@(x) median(x,2),PAT.FlyAmpMed.Vel,'UniformOutput',false))');
    PAT.GrandMed.Freq           = cell2mat((cellfun(@(x) median(x,2),PAT.FlyAmpMed.Freq,'UniformOutput',false))');
    PAT.GrandMed.Mag            = cell2mat((cellfun(@(x) median(x,2),PAT.FlyAmpMed.Mag,'UniformOutput',false))');
    PAT.GrandMed.Phase          = cell2mat((cellfun(@(x) circ_median(x,2),PAT.FlyAmpMed.Phase,'UniformOutput',false)))';

    WING.GrandMed.Time          = cell2mat((cellfun(@(x) median(x,2),WING.FlyAmpMed.Time,'UniformOutput',false))');
    WING.GrandMed.Pos           = cell2mat((cellfun(@(x) median(x,2),WING.FlyAmpMed.Pos,'UniformOutput',false))');
    WING.GrandMed.Vel           = cell2mat((cellfun(@(x) median(x,2),WING.FlyAmpMed.Vel,'UniformOutput',false))');
    WING.GrandMed.Freq          = cell2mat((cellfun(@(x) median(x,2),WING.FlyAmpMed.Freq,'UniformOutput',false))');
    WING.GrandMed.Mag           = cell2mat((cellfun(@(x) median(x,2),WING.FlyAmpMed.Mag,'UniformOutput',false))');
    WING.GrandMed.Phase      	= cell2mat((cellfun(@(x) circ_median(x,2),WING.FlyAmpMed.Phase,'UniformOutput',false)))';
    WING.GrandMed.COHR.Freq     = cell2mat((cellfun(@(x) median(x,2),WING.FlyAmpMed.COHR.Freq,'UniformOutput',false))');
    WING.GrandMed.COHR.Mag      = cell2mat((cellfun(@(x) median(x,2),WING.FlyAmpMed.COHR.Mag,'UniformOutput',false))');
    WING.GrandMed.GAIN      	= cell2mat((cellfun(@(x) median(x,2),WING.FlyAmpMed.GAIN,'UniformOutput',false))');
    WING.GrandMed.PHASE      	= cell2mat((cellfun(@(x) circ_median(x,2),WING.FlyAmpMed.PHASE,'UniformOutput',false)))';
    
% MEAN
%---------------------------------------------------------------------------------------------------------------------------------
    PAT.GrandMean.Time        	= cell2mat((cellfun(@(x) mean(x,2),PAT.FlyAmpMean.Time,'UniformOutput',false))');
    PAT.GrandMean.Pos         	= cell2mat((cellfun(@(x) mean(x,2),PAT.FlyAmpMean.Pos,'UniformOutput',false))');
    PAT.GrandMean.Vel        	= cell2mat((cellfun(@(x) mean(x,2),PAT.FlyAmpMean.Vel,'UniformOutput',false))');
    PAT.GrandMean.Freq       	= cell2mat((cellfun(@(x) mean(x,2),PAT.FlyAmpMean.Freq,'UniformOutput',false))');
    PAT.GrandMean.Mag         	= cell2mat((cellfun(@(x) mean(x,2),PAT.FlyAmpMean.Mag,'UniformOutput',false))');
    PAT.GrandMean.Phase       	= cell2mat((cellfun(@(x) circ_mean(x,[],2),PAT.FlyAmpMean.Phase,'UniformOutput',false)'));

    WING.GrandMean.Time         = cell2mat((cellfun(@(x) mean(x,2),WING.FlyAmpMean.Time,'UniformOutput',false))');
    WING.GrandMean.Pos      	= cell2mat((cellfun(@(x) mean(x,2),WING.FlyAmpMean.Pos,'UniformOutput',false))');
    WING.GrandMean.Vel          = cell2mat((cellfun(@(x) mean(x,2),WING.FlyAmpMean.Vel,'UniformOutput',false))');
    WING.GrandMean.Freq     	= cell2mat((cellfun(@(x) mean(x,2),WING.FlyAmpMean.Freq,'UniformOutput',false))');
    WING.GrandMean.Mag       	= cell2mat((cellfun(@(x) mean(x,2),WING.FlyAmpMean.Mag,'UniformOutput',false))');
    WING.GrandMean.Phase      	= cell2mat((cellfun(@(x) circ_mean(x,[],2),WING.FlyAmpMean.Phase,'UniformOutput',false)'));
    WING.GrandMean.COHR.Freq  	= cell2mat((cellfun(@(x) mean(x,2),WING.FlyAmpMean.COHR.Freq,'UniformOutput',false))');
    WING.GrandMean.COHR.Mag   	= cell2mat((cellfun(@(x) mean(x,2),WING.FlyAmpMean.COHR.Mag,'UniformOutput',false))');
    WING.GrandMean.GAIN      	= cell2mat((cellfun(@(x) mean(x,2),WING.FlyAmpMean.GAIN,'UniformOutput',false))');
    WING.GrandMean.PHASE      	= cell2mat((cellfun(@(x) circ_mean(x,[],2),WING.FlyAmpMean.PHASE,'UniformOutput',false)'));
    
% STD
%---------------------------------------------------------------------------------------------------------------------------------
    PAT.GrandSTD.Time           = cell2mat((cellfun(@(x) std(x,0,2),PAT.FlyAmpSTD.Time,'UniformOutput',false))');
    PAT.GrandSTD.Pos            = cell2mat((cellfun(@(x) std(x,0,2),PAT.FlyAmpSTD.Pos,'UniformOutput',false))');
    PAT.GrandSTD.Vel            = cell2mat((cellfun(@(x) std(x,0,2),PAT.FlyAmpSTD.Vel,'UniformOutput',false))');
    PAT.GrandSTD.Freq           = cell2mat((cellfun(@(x) std(x,0,2),PAT.FlyAmpSTD.Freq,'UniformOutput',false))');
    PAT.GrandSTD.Mag            = cell2mat((cellfun(@(x) std(x,0,2),PAT.FlyAmpSTD.Mag,'UniformOutput',false))');
    PAT.GrandSTD.Phase          = cell2mat((cellfun(@(x) circ_std(x,[],[],2),PAT.FlyAmpSTD.Phase,'UniformOutput',false))');

    WING.GrandSTD.Time          = cell2mat((cellfun(@(x) std(x,0,2),WING.FlyAmpSTD.Time,'UniformOutput',false))');
    WING.GrandSTD.Pos           = cell2mat((cellfun(@(x) std(x,0,2),WING.FlyAmpSTD.Pos,'UniformOutput',false))');
    WING.GrandSTD.Vel           = cell2mat((cellfun(@(x) std(x,0,2),WING.FlyAmpSTD.Vel,'UniformOutput',false))');
    WING.GrandSTD.Freq          = cell2mat((cellfun(@(x) std(x,0,2),WING.FlyAmpSTD.Freq,'UniformOutput',false))');
    WING.GrandSTD.Mag           = cell2mat((cellfun(@(x) std(x,0,2),WING.FlyAmpSTD.Mag,'UniformOutput',false))');
    WING.GrandSTD.Phase      	= cell2mat((cellfun(@(x) circ_std(x,[],[],2),WING.FlyAmpSTD.Phase,'UniformOutput',false))');
    WING.GrandSTD.COHR.Freq     = cell2mat((cellfun(@(x) std(x,0,2),WING.FlyAmpSTD.COHR.Freq,'UniformOutput',false))');
    WING.GrandSTD.COHR.Mag      = cell2mat((cellfun(@(x) std(x,0,2),WING.FlyAmpSTD.COHR.Mag,'UniformOutput',false))');
    WING.GrandSTD.GAIN      	= cell2mat((cellfun(@(x) std(x,0,2),WING.FlyAmpSTD.GAIN,'UniformOutput',false))');
    WING.GrandSTD.PHASE      	= cell2mat((cellfun(@(x) circ_std(x,[],[],2),WING.FlyAmpSTD.PHASE,'UniformOutput',false))');
    
%% Save ouputs as structure %%
%---------------------------------------------------------------------------------------------------------------------------------
disp('Saving...')
save([rootdir 'DATA\' filename '.mat'],'PAT','WING','FD','T','n','unq')
disp('DONE')
end