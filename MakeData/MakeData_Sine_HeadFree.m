function [PAT,WING,HEAD,BODE,FD,T,n,unq] = MakeData_Sine_HeadFree(rootdir,Amp)
%% MakeData_Sine_HeadFree: Reads in all raw trials, transforms data, and saves in organized structure for use with figure functions
%   INPUTS:
%       root        : root directory
%       Amp         : sinusoid amplitude, set's directory
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
% Amp = 3.75;
% rootdir = 'H:\EXPERIMENTS\Experiment_Sinusoid\';
%---------------------------------------------------------------------------------------------------------------------------------
%% Setup Directories %%
%---------------------------------------------------------------------------------------------------------------------------------
filename = ['Sine_HeadFree_' num2str(Amp) '_DATA_' ] ;
root.daq = ['H:\EXPERIMENTS\Experiment_Sinusoid\' num2str(Amp) '\'];
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
FD.Freq      = zeros(n.Trial,1);     % amplitude [deg]
for jj = 1:n.Trial
    temp = textscan(char(FILES{jj}), '%s', 'delimiter', '_.'); temp = temp{1} ; % read individual strings into temp variable
    FD.Fly(jj,1)    = str2double(temp{2}); % store fly #
    FD.Trial(jj,1)  = str2double(temp{4}); % store trial #
    
    if ~isnan(str2double(temp{7}))
        FD.Freq(jj,1) = str2double([temp{6} '.' temp{7}]); % store amplitude
    else
        FD.Freq(jj,1) = str2double(temp{6});
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
unq.Freq 	= sort(unique(FD.Freq));	% find all unique amplitudes
n.Freq       = length(unq.Freq);  	% # of unique amplitudes
FD.idxFreq   = zeros(n.Trial,1);     % amplitude index
for kk = 1:n.Freq
   FD.idxFreq(FD.Freq == unq.Freq(kk)) = kk; % frequency index
end

fprintf('Total Flies''s: %i \n',n.Fly) ; fprintf('Total Trials''s: %i \n',n.Trial)
T = cell2table(FD.trialFly,'VariableNames',{'Fly','Trials'});
disp(T)
clear pp kk
%% Get Data %%
%---------------------------------------------------------------------------------------------------------------------------------
disp('Loading...')
% Preallocate data cells
WING.ALL.Pos = cell(n.Freq,1);
HEAD.ALL.Pos = cell(n.Freq,1);
for kk = 1:n.Fly
    for jj = 1:n.Freq
        % PATTERN Data
        PAT.Pos{kk,1}{jj,1}         = [];
     	PAT.Time{kk,1}{jj,1}        = [];
        PAT.Vel{kk,1}{jj,1}         = [];
        PAT.VelMean{kk,1}{jj,1} 	= [];
        PAT.VelSTD{kk,1}{jj,1}      = [];
        PAT.Freq{kk,1}{jj,1}        = [];
        PAT.Mag{kk,1}{jj,1}     	= [];
        PAT.Phase{kk,1}{jj,1}       = [];
        PAT.MAG{kk,1}{jj,1}     	= [];
        PAT.PHASE{kk,1}{jj,1}       = [];
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
        WING.MAG{kk,1}{jj,1}        = [];
        WING.PHASE{kk,1}{jj,1}      = [];
        WING.GAIN{kk,1}{jj,1}       = [];
        WING.PhaseDiff{kk,1}{jj,1} 	= [];
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
        HEAD.Err.MAG{kk,1}{jj,1}	= [];
        HEAD.Err.PHASE{kk,1}{jj,1}	= [];
    	HEAD.COHR.Freq{kk,1}{jj,1}  = [];
        HEAD.COHR.Mag{kk,1}{jj,1}   = [];
        HEAD.MAG{kk,1}{jj,1}        = [];
        HEAD.PHASE{kk,1}{jj,1}      = [];
        HEAD.GAIN{kk,1}{jj,1}       = [];
        HEAD.PhaseDiff{kk,1}{jj,1} 	= [];
        % BODE
        BODE.head2wing.GAIN{kk,1}{jj,1}         = [];
        BODE.head2wing.PhaseDiff{kk,1}{jj,1}	= [];
    end
end
% Store data in organized cells
tt = (0:1/200:10)';
tt = tt(1:end-1);
count = 0;
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
        count = count + 1;
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
    wing.Pos    = interp1(wing.Time, wing.Pos, tt, 'nearest'); % interpolate pattern to match fly video [V]
    wing.Vel    = interp1(wing.Time, wing.Vel, tt, 'nearest'); % interpolate pattern to match fly video [V/s]
    wing.Time  	= tt; % new wing time [s]
    wing.Fs  	= head.Fs ; % new wing sampling frequency [Hz]
	%-----------------------------------------------------------------------------------------------------------------------------
 	head.Err.Pos    = pat.Pos - head.Pos; % calculate Error between head & pattern (retinal slip) [deg]
    head.Err.Vel    = pat.Vel - head.Vel; % calculate Error between head & pattern (retinal slip) [deg/s]
    %-----------------------------------------------------------------------------------------------------------------------------
    % Convert head, wings, & pattern data to frequency domain using FFT %
    [head.Freq , head.Mag , head.Phase]	= FFT(head.Time,head.Pos);
    [wing.Freq , wing.Mag , wing.Phase]	= FFT(wing.Time,wing.Pos);
    [pat.Freq  , pat.Mag  , pat.Phase ] = FFT(pat.Time,pat.Pos);
	[head.Err.Freq , head.Err.Mag, head.Err.Phase]  = FFT(head.Time,head.Err.Pos);
	%-----------------------------------------------------------------------------------------------------------------------------
    % Calculate magnitude & phase at input frequency %
    [head.MAG,head.PHASE] = Get_IO_Freq(head.Freq,head.Mag,head.Phase,FD.Freq(kk),0.2,false);
    [wing.MAG,wing.PHASE] = Get_IO_Freq(wing.Freq,wing.Mag,wing.Phase,FD.Freq(kk));
 	[pat.MAG,pat.PHASE] = Get_IO_Freq(pat.Freq,pat.Mag,pat.Phase,FD.Freq(kk),0.2,false);
  	[head.Err.MAG,head.Err.PHASE] = Get_IO_Freq(head.Err.Freq,head.Err.Mag,head.Err.Phase,FD.Freq(kk));
	%-----------------------------------------------------------------------------------------------------------------------------
    % Calculate BODE at input frequency %
    head.GAIN       = head.MAG./pat.MAG;
    head.PhaseDiff  = head.PHASE - pat.PHASE;
	wing.GAIN       = wing.MAG./head.Err.MAG;
    wing.PhaseDiff  = wing.PHASE - head.Err.PHASE;
    bode.head2wing.Gain         = wing.MAG./head.MAG;
    bode.head2wing.PhaseDiff    = wing.PHASE - head.PHASE;
    %-----------------------------------------------------------------------------------------------------------------------------
    % Calculate coherence %
    [head.cohr.mag,head.cohr.f] = mscohere(pat.Pos , head.Pos ,[],[] , head.Freq , head.Fs);
	[wing.cohr.mag,wing.cohr.f] = mscohere(pat.Pos , wing.Pos ,[],[] , wing.Freq , wing.Fs);
    %-----------------------------------------------------------------------------------------------------------------------------
    % Store data in cells %
  	% PATTERN
	PAT.Time        {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = pat.Time;
	PAT.Pos         {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = pat.Pos;
	PAT.Vel         {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = pat.Vel;
   	PAT.VelMean     {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = pat.VelMean;
	PAT.VelSTD      {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = pat.VelSTD;
    PAT.Freq        {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = pat.Freq;
	PAT.Mag         {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = pat.Mag;
	PAT.Phase       {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = pat.Phase;
	PAT.MAG         {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = pat.MAG;
	PAT.PHASE       {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = pat.PHASE;
    % WINGS
    WING.ALL.Pos    {FD.idxFreq(kk)}(:,end+1) = wing.Pos;
	WING.Pos        {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = wing.Pos;
	WING.Time       {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = wing.Time;
	WING.Vel        {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = wing.Vel;
	WING.VelMean	{FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = wing.VelMean;
	WING.VelSTD     {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = wing.VelSTD;
	WING.Freq       {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = wing.Freq;
	WING.Mag        {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = wing.Mag;
	WING.Phase      {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = wing.Phase;
	WING.MAG      	{FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = wing.MAG;
	WING.PHASE    	{FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = wing.PHASE;
	WING.COHR.Freq  {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = wing.cohr.f;
    WING.COHR.Mag   {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = wing.cohr.mag;
	WING.GAIN      	{FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = wing.GAIN;
	WING.PhaseDiff	{FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = wing.PhaseDiff;
    % HEAD
    HEAD.ALL.Pos    {FD.idxFreq(kk)}(:,end+1) = head.Pos;
	HEAD.Time       {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = head.Time;
	HEAD.Pos        {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = head.Pos;
    HEAD.Vel        {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = head.Vel;
    HEAD.VelMean    {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = head.VelMean;
    HEAD.VelSTD     {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = head.VelSTD;
	HEAD.Err.Pos    {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = head.Err.Pos;
    HEAD.Err.Freq 	{FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = head.Err.Freq;
    HEAD.Err.Mag    {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = head.Err.Mag;
 	HEAD.Err.Phase 	{FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = head.Err.Phase;
	HEAD.Freq       {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = head.Freq;
	HEAD.Mag        {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = head.Mag;
	HEAD.Phase      {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = head.Phase;
	HEAD.MAG      	{FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = head.MAG;
	HEAD.PHASE    	{FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = head.PHASE;
	HEAD.COHR.Freq  {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = head.cohr.f;
    HEAD.COHR.Mag   {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = head.cohr.mag;
	HEAD.GAIN      	{FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = head.GAIN;
	HEAD.PhaseDiff	{FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = head.PhaseDiff;
    % BODE
    BODE.head2wing.GAIN         {FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = bode.head2wing.Gain;
	BODE.head2wing.PhaseDiff  	{FD.idxFly(kk),1}{FD.idxFreq(kk),1}(:,end+1) = bode.head2wing.PhaseDiff;
end
clear jj kk a b t_p t_v hAngles data  wing pat bode tt
disp('LOADING DONE')
disp('Bad Trial:')
disp(count)
pause(3)
%% FLY Stats by Fly %%
%---------------------------------------------------------------------------------------------------------------------------------
for kk = 1:n.Fly
    for jj = 1:n.Freq
    % MEDIAN
    %-----------------------------------------------------------------------------------------------------------------------------
        % PATTERN
        PAT.FlyMed.Time         {kk,1}(:,jj)	= median(PAT.Time{kk}{jj},2);
        PAT.FlyMed.Pos          {kk,1}(:,jj)	= median(PAT.Pos{kk}{jj},2);
        PAT.FlyMed.Vel          {kk,1}(:,jj)	= median(PAT.Vel{kk}{jj},2);
        PAT.FlyMed.Freq         {kk,1}(:,jj)	= median(PAT.Freq{kk}{jj},2);
        PAT.FlyMed.Mag          {kk,1}(:,jj)	= median(PAT.Mag{kk}{jj},2);
        PAT.FlyMed.Phase     	{kk,1}(:,jj)	= circ_median(PAT.Phase{kk}{jj},2)';
        % HEAD
        HEAD.FlyMed.Time      	{kk,1}(:,jj)	= median(HEAD.Time{kk}{jj},2);
        HEAD.FlyMed.Pos         {kk,1}(:,jj)	= median(HEAD.Pos{kk}{jj},2);
        HEAD.FlyMed.Vel         {kk,1}(:,jj)	= median(HEAD.Vel{kk}{jj},2);
        HEAD.FlyMed.Freq        {kk,1}(:,jj)	= median(HEAD.Freq{kk}{jj},2);
        HEAD.FlyMed.Mag         {kk,1}(:,jj)	= median(HEAD.Mag{kk}{jj},2);
        HEAD.FlyMed.Phase       {kk,1}(:,jj)	= circ_median(HEAD.Phase{kk}{jj},2)';
        HEAD.FlyMed.COHR.Freq 	{kk,1}(:,jj)	= median(HEAD.COHR.Freq{kk}{jj},2);
        HEAD.FlyMed.COHR.Mag 	{kk,1}(:,jj)	= median(HEAD.COHR.Mag{kk}{jj},2);
        HEAD.FlyMed.GAIN      	{kk,1}(:,jj)	= median(HEAD.GAIN{kk}{jj},2);
        HEAD.FlyMed.PhaseDiff  	{kk,1}(:,jj)	= circ_median(HEAD.PhaseDiff{kk}{jj},2)';

        HEAD.FlyMed.Err.Pos    	{kk,1}(:,jj)	= median(HEAD.Err.Pos{kk}{jj},2);
        HEAD.FlyMed.Err.Freq  	{kk,1}(:,jj)	= median(HEAD.Err.Freq{kk}{jj},2);
        HEAD.FlyMed.Err.Mag    	{kk,1}(:,jj)	= median(HEAD.Err.Mag{kk}{jj},2);
    	HEAD.FlyMed.Err.Phase 	{kk,1}(:,jj)	= circ_median(HEAD.Err.Phase{kk}{jj},2)';
        % WINGS
        WING.FlyMed.Time    	{kk,1}(:,jj)	= median(WING.Time{kk}{jj},2);
        WING.FlyMed.Pos         {kk,1}(:,jj)	= median(WING.Pos{kk}{jj},2);
        WING.FlyMed.Vel         {kk,1}(:,jj)	= median(WING.Vel{kk}{jj},2);
        WING.FlyMed.Freq    	{kk,1}(:,jj)	= median(WING.Freq{kk}{jj},2);
        WING.FlyMed.Mag     	{kk,1}(:,jj)	= median(WING.Mag{kk}{jj},2);
        WING.FlyMed.Phase    	{kk,1}(:,jj)	= circ_median(WING.Phase{kk}{jj},2)';
        WING.FlyMed.COHR.Freq 	{kk,1}(:,jj)	= median(WING.COHR.Freq{kk}{jj},2);
        WING.FlyMed.COHR.Mag 	{kk,1}(:,jj)	= median(WING.COHR.Mag{kk}{jj},2);
        WING.FlyMed.GAIN      	{kk,1}(:,jj)	= median(WING.GAIN{kk}{jj},2);
        WING.FlyMed.PhaseDiff  	{kk,1}(:,jj)	= circ_median(WING.PhaseDiff{kk}{jj},2)';
        % BODE
     	BODE.FlyMed.head2wing.GAIN          {kk,1}(:,jj)	= median(BODE.head2wing.GAIN{kk}{jj},2);
     	BODE.FlyMed.head2wing.PhaseDiff     {kk,1}(:,jj) 	= circ_median(BODE.head2wing.PhaseDiff{kk}{jj},2)';
	% MEAN
    %-----------------------------------------------------------------------------------------------------------------------------
        % PATTERN
        PAT.FlyMean.Time    	{kk,1}(:,jj)	= mean(PAT.Time{kk}{jj},2);
        PAT.FlyMean.Pos       	{kk,1}(:,jj)	= mean(PAT.Pos{kk}{jj},2);
        PAT.FlyMean.Vel      	{kk,1}(:,jj)	= mean(PAT.Vel{kk}{jj},2);
        PAT.FlyMean.Freq     	{kk,1}(:,jj)	= mean(PAT.Freq{kk}{jj},2);
        PAT.FlyMean.Mag       	{kk,1}(:,jj)	= mean(PAT.Mag{kk}{jj},2);
        PAT.FlyMean.Phase     	{kk,1}(:,jj)	= circ_mean(PAT.Phase{kk}{jj},[],2)';
        % HEAD
        HEAD.FlyMean.Time      	{kk,1}(:,jj)	= mean(HEAD.Time{kk}{jj},2);
        HEAD.FlyMean.Pos      	{kk,1}(:,jj)	= mean(HEAD.Pos{kk}{jj},2);
        HEAD.FlyMean.Vel     	{kk,1}(:,jj)	= mean(HEAD.Vel{kk}{jj},2);
        HEAD.FlyMean.Freq     	{kk,1}(:,jj)	= mean(HEAD.Freq{kk}{jj},2);
        HEAD.FlyMean.Mag    	{kk,1}(:,jj)	= mean(HEAD.Mag{kk}{jj},2);
        HEAD.FlyMean.Phase    	{kk,1}(:,jj)	= circ_mean(HEAD.Phase{kk}{jj},[],2)';
        HEAD.FlyMean.COHR.Freq 	{kk,1}(:,jj)	= mean(HEAD.COHR.Freq{kk}{jj},2);
        HEAD.FlyMean.COHR.Mag 	{kk,1}(:,jj)	= mean(HEAD.COHR.Mag{kk}{jj},2);
        HEAD.FlyMean.GAIN   	{kk,1}(:,jj)	= mean(HEAD.GAIN{kk}{jj},2);
        HEAD.FlyMean.PhaseDiff 	{kk,1}(:,jj)	= circ_mean(HEAD.PhaseDiff{kk}{jj},[],2)';

        HEAD.FlyMean.Err.Pos 	{kk,1}(:,jj)	= mean(HEAD.Err.Pos{kk}{jj},2);
        HEAD.FlyMean.Err.Freq  	{kk,1}(:,jj)	= mean(HEAD.Err.Freq{kk}{jj},2);
        HEAD.FlyMean.Err.Mag   	{kk,1}(:,jj)	= mean(HEAD.Err.Mag{kk}{jj},2);
    	HEAD.FlyMean.Err.Phase 	{kk,1}(:,jj)	= circ_mean(HEAD.Err.Phase{kk}{jj},[],2)';
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
        WING.FlyMean.PhaseDiff 	{kk,1}(:,jj)	= circ_mean(WING.PhaseDiff{kk}{jj},[],2)';
        % BODE
     	BODE.FlyMean.head2wing.GAIN         {kk,1}(:,jj)	= mean(BODE.head2wing.GAIN{kk}{jj},2);
     	BODE.FlyMean.head2wing.PhaseDiff	{kk,1}(:,jj) 	= circ_mean(BODE.head2wing.PhaseDiff{kk}{jj},[],2)';        
    % STD
    %-----------------------------------------------------------------------------------------------------------------------------
        % PATTERN
        PAT.FlySTD.Time         {kk,1}(:,jj)	= std(PAT.Time{kk}{jj},0,2);
        PAT.FlySTD.Pos          {kk,1}(:,jj)	= std(PAT.Pos{kk}{jj},0,2);
        PAT.FlySTD.Vel          {kk,1}(:,jj)	= std(PAT.Vel{kk}{jj},0,2);
        PAT.FlySTD.Freq         {kk,1}(:,jj)	= std(PAT.Freq{kk}{jj},0,2);
        PAT.FlySTD.Mag          {kk,1}(:,jj)	= std(PAT.Mag{kk}{jj},0,2);
        PAT.FlySTD.Phase     	{kk,1}(:,jj)	= circ_std(PAT.Phase{kk}{jj},[],[],2);
        % HEAD
        HEAD.FlySTD.Time      	{kk,1}(:,jj)	= std(HEAD.Time{kk}{jj},0,2);
        HEAD.FlySTD.Pos         {kk,1}(:,jj)	= std(HEAD.Pos{kk}{jj},0,2);
        HEAD.FlySTD.Vel         {kk,1}(:,jj)	= std(HEAD.Vel{kk}{jj},0,2);
        HEAD.FlySTD.Freq        {kk,1}(:,jj)	= std(HEAD.Freq{kk}{jj},0,2);
        HEAD.FlySTD.Mag         {kk,1}(:,jj)	= std(HEAD.Mag{kk}{jj},0,2);
        HEAD.FlySTD.Phase       {kk,1}(:,jj)	= circ_std(HEAD.Phase{kk}{jj},[],[],2);
        HEAD.FlySTD.COHR.Freq 	{kk,1}(:,jj)	= std(HEAD.COHR.Freq{kk}{jj},0,2);
        HEAD.FlySTD.COHR.Mag 	{kk,1}(:,jj)	= std(HEAD.COHR.Mag{kk}{jj},0,2);
        HEAD.FlySTD.GAIN        {kk,1}(:,jj)	= std(HEAD.GAIN{kk}{jj},0,2);
        HEAD.FlySTD.PhaseDiff  	{kk,1}(:,jj)	= circ_std(HEAD.PhaseDiff{kk}{jj},[],[],2);

        HEAD.FlySTD.Err.Pos    	{kk,1}(:,jj)	= std(HEAD.Err.Pos{kk}{jj},0,2);
        HEAD.FlySTD.Err.Freq  	{kk,1}(:,jj)	= std(HEAD.Err.Freq{kk}{jj},0,2);
        HEAD.FlySTD.Err.Mag    	{kk,1}(:,jj)	= std(HEAD.Err.Mag{kk}{jj},0,2);
    	HEAD.FlySTD.Err.Phase 	{kk,1}(:,jj)	= circ_std(HEAD.Err.Phase{kk}{jj},[],[],2);
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
        WING.FlySTD.PhaseDiff  	{kk,1}(:,jj)	= circ_std(WING.PhaseDiff{kk}{jj},[],[],2);
        % BODE
     	BODE.FlySTD.head2wing.GAIN          {kk,1}(:,jj)	= std(BODE.head2wing.GAIN{kk}{jj},0,2);
     	BODE.FlySTD.head2wing.PhaseDiff     {kk,1}(:,jj) 	= circ_std(BODE.head2wing.PhaseDiff{kk}{jj},[],[],2);           
    end
end
clear kk jj
%% FLY stats by Amplitude %%
%---------------------------------------------------------------------------------------------------------------------------------
for jj = 1:n.Freq
    for kk = 1:n.Fly
	% MEDIAN
    %-----------------------------------------------------------------------------------------------------------------------------
        PAT.FlyAmpMed.Time          {jj,1}(:,kk)  	= PAT.FlyMed.Time{kk}(:,jj);
        PAT.FlyAmpMed.Pos           {jj,1}(:,kk) 	= PAT.FlyMed.Pos{kk}(:,jj);
    	PAT.FlyAmpMed.Vel           {jj,1}(:,kk)  	= PAT.FlyMed.Vel{kk}(:,jj);
        PAT.FlyAmpMed.Freq          {jj,1}(:,kk)  	= PAT.FlyMed.Freq{kk}(:,jj);
        PAT.FlyAmpMed.Phase         {jj,1}(:,kk)  	= PAT.FlyMed.Phase{kk}(:,jj);
        PAT.FlyAmpMed.Mag           {jj,1}(:,kk)   	= PAT.FlyMed.Mag{kk}(:,jj);
        
        HEAD.FlyAmpMed.Time         {jj,1}(:,kk)  	= HEAD.FlyMed.Time{kk}(:,jj);
        HEAD.FlyAmpMed.Pos          {jj,1}(:,kk)   	= HEAD.FlyMed.Pos{kk}(:,jj);
    	HEAD.FlyAmpMed.Vel          {jj,1}(:,kk)   	= HEAD.FlyMed.Vel{kk}(:,jj);
        HEAD.FlyAmpMed.Freq         {jj,1}(:,kk)   	= HEAD.FlyMed.Freq{kk}(:,jj);
        HEAD.FlyAmpMed.Phase        {jj,1}(:,kk)  	= HEAD.FlyMed.Phase{kk}(:,jj);
        HEAD.FlyAmpMed.Mag          {jj,1}(:,kk)   	= HEAD.FlyMed.Mag{kk}(:,jj);
        HEAD.FlyAmpMed.COHR.Freq    {jj,1}(:,kk)  	= HEAD.FlyMed.COHR.Freq{kk}(:,jj);
        HEAD.FlyAmpMed.COHR.Mag     {jj,1}(:,kk) 	= HEAD.FlyMed.COHR.Mag{kk}(:,jj);
        HEAD.FlyAmpMed.GAIN         {jj,1}(:,kk)  	= HEAD.FlyMed.GAIN{kk}(:,jj);
        HEAD.FlyAmpMed.PhaseDiff  	{jj,1}(:,kk)   	= HEAD.FlyMed.PhaseDiff{kk}(:,jj);
        HEAD.FlyAmpMed.Err.Pos      {jj,1}(:,kk) 	= HEAD.FlyMed.Err.Pos{kk}(:,jj);
    	HEAD.FlyAmpMed.Err.Freq     {jj,1}(:,kk) 	= HEAD.FlyMed.Err.Freq{kk}(:,jj);
    	HEAD.FlyAmpMed.Err.Mag      {jj,1}(:,kk) 	= HEAD.FlyMed.Err.Mag{kk}(:,jj);
    	HEAD.FlyAmpMed.Err.Phase    {jj,1}(:,kk) 	= HEAD.FlyMed.Err.Phase{kk}(:,jj);
        
        WING.FlyAmpMed.Time         {jj,1}(:,kk)  	= WING.FlyMed.Time{kk}(:,jj);
        WING.FlyAmpMed.Pos          {jj,1}(:,kk)   	= WING.FlyMed.Pos{kk}(:,jj);
    	WING.FlyAmpMed.Vel          {jj,1}(:,kk)   	= WING.FlyMed.Vel{kk}(:,jj);
        WING.FlyAmpMed.Freq         {jj,1}(:,kk)   	= WING.FlyMed.Freq{kk}(:,jj);
        WING.FlyAmpMed.Phase        {jj,1}(:,kk)  	= WING.FlyMed.Phase{kk}(:,jj);
        WING.FlyAmpMed.Mag          {jj,1}(:,kk)  	= WING.FlyMed.Mag{kk}(:,jj);
        WING.FlyAmpMed.COHR.Freq    {jj,1}(:,kk)	= WING.FlyMed.COHR.Freq{kk}(:,jj);
        WING.FlyAmpMed.COHR.Mag     {jj,1}(:,kk)	= WING.FlyMed.COHR.Mag{kk}(:,jj);
        WING.FlyAmpMed.GAIN         {jj,1}(:,kk)   	= WING.FlyMed.GAIN{kk}(:,jj);
        WING.FlyAmpMed.PhaseDiff  	{jj,1}(:,kk)   	= WING.FlyMed.PhaseDiff{kk}(:,jj);
        
        BODE.FlyAmpMed.head2wing.GAIN       {jj,1}(:,kk)    = BODE.FlyMed.head2wing.GAIN{kk}(:,jj);
        BODE.FlyAmpMed.head2wing.PhaseDiff  {jj,1}(:,kk)    = BODE.FlyMed.head2wing.PhaseDiff{kk}(:,jj);
	% MEAN
    %-----------------------------------------------------------------------------------------------------------------------------
        PAT.FlyAmpMean.Time        	{jj,1}(:,kk)  	= PAT.FlyMean.Time{kk}(:,jj);
        PAT.FlyAmpMean.Pos         	{jj,1}(:,kk) 	= PAT.FlyMean.Pos{kk}(:,jj);
    	PAT.FlyAmpMean.Vel         	{jj,1}(:,kk)  	= PAT.FlyMean.Vel{kk}(:,jj);
        PAT.FlyAmpMean.Freq      	{jj,1}(:,kk)  	= PAT.FlyMean.Freq{kk}(:,jj);
        PAT.FlyAmpMean.Phase     	{jj,1}(:,kk)  	= PAT.FlyMean.Phase{kk}(:,jj);
        PAT.FlyAmpMean.Mag       	{jj,1}(:,kk)   	= PAT.FlyMean.Mag{kk}(:,jj);
        
        HEAD.FlyAmpMean.Time     	{jj,1}(:,kk)  	= HEAD.FlyMean.Time{kk}(:,jj);
        HEAD.FlyAmpMean.Pos      	{jj,1}(:,kk)   	= HEAD.FlyMean.Pos{kk}(:,jj);
    	HEAD.FlyAmpMean.Vel     	{jj,1}(:,kk)   	= HEAD.FlyMean.Vel{kk}(:,jj);
        HEAD.FlyAmpMean.Freq      	{jj,1}(:,kk)   	= HEAD.FlyMean.Freq{kk}(:,jj);
        HEAD.FlyAmpMean.Phase      	{jj,1}(:,kk)  	= HEAD.FlyMean.Phase{kk}(:,jj);
        HEAD.FlyAmpMean.Mag      	{jj,1}(:,kk)   	= HEAD.FlyMean.Mag{kk}(:,jj);
        HEAD.FlyAmpMean.COHR.Freq  	{jj,1}(:,kk)  	= HEAD.FlyMean.COHR.Freq{kk}(:,jj);
        HEAD.FlyAmpMean.COHR.Mag   	{jj,1}(:,kk) 	= HEAD.FlyMean.COHR.Mag{kk}(:,jj);
        HEAD.FlyAmpMean.GAIN       	{jj,1}(:,kk)  	= HEAD.FlyMean.GAIN{kk}(:,jj);
        HEAD.FlyAmpMean.PhaseDiff 	{jj,1}(:,kk)   	= HEAD.FlyMean.PhaseDiff{kk}(:,jj);
        HEAD.FlyAmpMean.Err.Pos   	{jj,1}(:,kk) 	= HEAD.FlyMean.Err.Pos{kk}(:,jj);
    	HEAD.FlyAmpMean.Err.Freq  	{jj,1}(:,kk) 	= HEAD.FlyMean.Err.Freq{kk}(:,jj);
    	HEAD.FlyAmpMean.Err.Mag 	{jj,1}(:,kk) 	= HEAD.FlyMean.Err.Mag{kk}(:,jj);
    	HEAD.FlyAmpMean.Err.Phase  	{jj,1}(:,kk) 	= HEAD.FlyMean.Err.Phase{kk}(:,jj);
        
        WING.FlyAmpMean.Time     	{jj,1}(:,kk)  	= WING.FlyMean.Time{kk}(:,jj);
        WING.FlyAmpMean.Pos      	{jj,1}(:,kk)   	= WING.FlyMean.Pos{kk}(:,jj);
    	WING.FlyAmpMean.Vel        	{jj,1}(:,kk)   	= WING.FlyMean.Vel{kk}(:,jj);
        WING.FlyAmpMean.Freq    	{jj,1}(:,kk)   	= WING.FlyMean.Freq{kk}(:,jj);
        WING.FlyAmpMean.Phase    	{jj,1}(:,kk)  	= WING.FlyMean.Phase{kk}(:,jj);
        WING.FlyAmpMean.Mag       	{jj,1}(:,kk)  	= WING.FlyMean.Mag{kk}(:,jj);
        WING.FlyAmpMean.COHR.Freq  	{jj,1}(:,kk)	= WING.FlyMean.COHR.Freq{kk}(:,jj);
        WING.FlyAmpMean.COHR.Mag  	{jj,1}(:,kk)	= WING.FlyMean.COHR.Mag{kk}(:,jj);
        WING.FlyAmpMean.GAIN      	{jj,1}(:,kk)   	= WING.FlyMean.GAIN{kk}(:,jj);
        WING.FlyAmpMean.PhaseDiff  	{jj,1}(:,kk)   	= WING.FlyMean.PhaseDiff{kk}(:,jj);
        
        BODE.FlyAmpMean.head2wing.GAIN          {jj,1}(:,kk)    = BODE.FlyMean.head2wing.GAIN{kk}(:,jj);
        BODE.FlyAmpMean.head2wing.PhaseDiff     {jj,1}(:,kk)    = BODE.FlyMean.head2wing.PhaseDiff{kk}(:,jj);
	% STD
    %-----------------------------------------------------------------------------------------------------------------------------
        PAT.FlyAmpSTD.Time        	{jj,1}(:,kk)  	= PAT.FlySTD.Time{kk}(:,jj);
        PAT.FlyAmpSTD.Pos         	{jj,1}(:,kk) 	= PAT.FlySTD.Pos{kk}(:,jj);
    	PAT.FlyAmpSTD.Vel         	{jj,1}(:,kk)  	= PAT.FlySTD.Vel{kk}(:,jj);
        PAT.FlyAmpSTD.Freq      	{jj,1}(:,kk)  	= PAT.FlySTD.Freq{kk}(:,jj);
        PAT.FlyAmpSTD.Phase     	{jj,1}(:,kk)  	= PAT.FlySTD.Phase{kk}(:,jj);
        PAT.FlyAmpSTD.Mag       	{jj,1}(:,kk)   	= PAT.FlySTD.Mag{kk}(:,jj);
        
        HEAD.FlyAmpSTD.Time     	{jj,1}(:,kk)  	= HEAD.FlySTD.Time{kk}(:,jj);
        HEAD.FlyAmpSTD.Pos      	{jj,1}(:,kk)   	= HEAD.FlySTD.Pos{kk}(:,jj);
    	HEAD.FlyAmpSTD.Vel          {jj,1}(:,kk)   	= HEAD.FlySTD.Vel{kk}(:,jj);
        HEAD.FlyAmpSTD.Freq      	{jj,1}(:,kk)   	= HEAD.FlySTD.Freq{kk}(:,jj);
        HEAD.FlyAmpSTD.Phase      	{jj,1}(:,kk)  	= HEAD.FlySTD.Phase{kk}(:,jj);
        HEAD.FlyAmpSTD.Mag      	{jj,1}(:,kk)   	= HEAD.FlySTD.Mag{kk}(:,jj);
        HEAD.FlyAmpSTD.COHR.Freq  	{jj,1}(:,kk)  	= HEAD.FlySTD.COHR.Freq{kk}(:,jj);
        HEAD.FlyAmpSTD.COHR.Mag   	{jj,1}(:,kk) 	= HEAD.FlySTD.COHR.Mag{kk}(:,jj);
        HEAD.FlyAmpSTD.GAIN       	{jj,1}(:,kk)  	= HEAD.FlySTD.GAIN{kk}(:,jj);
        HEAD.FlyAmpSTD.PhaseDiff    {jj,1}(:,kk)   	= HEAD.FlySTD.PhaseDiff{kk}(:,jj);
        HEAD.FlyAmpSTD.Err.Pos   	{jj,1}(:,kk) 	= HEAD.FlySTD.Err.Pos{kk}(:,jj);
    	HEAD.FlyAmpSTD.Err.Freq  	{jj,1}(:,kk) 	= HEAD.FlySTD.Err.Freq{kk}(:,jj);
    	HEAD.FlyAmpSTD.Err.Mag      {jj,1}(:,kk) 	= HEAD.FlySTD.Err.Mag{kk}(:,jj);
    	HEAD.FlyAmpSTD.Err.Phase  	{jj,1}(:,kk) 	= HEAD.FlySTD.Err.Phase{kk}(:,jj);
        
        WING.FlyAmpSTD.Time     	{jj,1}(:,kk)  	= WING.FlySTD.Time{kk}(:,jj);
        WING.FlyAmpSTD.Pos      	{jj,1}(:,kk)   	= WING.FlySTD.Pos{kk}(:,jj);
    	WING.FlyAmpSTD.Vel        	{jj,1}(:,kk)   	= WING.FlySTD.Vel{kk}(:,jj);
        WING.FlyAmpSTD.Freq         {jj,1}(:,kk)   	= WING.FlySTD.Freq{kk}(:,jj);
        WING.FlyAmpSTD.Phase    	{jj,1}(:,kk)  	= WING.FlySTD.Phase{kk}(:,jj);
        WING.FlyAmpSTD.Mag       	{jj,1}(:,kk)  	= WING.FlySTD.Mag{kk}(:,jj);
        WING.FlyAmpSTD.COHR.Freq  	{jj,1}(:,kk)	= WING.FlySTD.COHR.Freq{kk}(:,jj);
        WING.FlyAmpSTD.COHR.Mag  	{jj,1}(:,kk)	= WING.FlySTD.COHR.Mag{kk}(:,jj);
        WING.FlyAmpSTD.GAIN      	{jj,1}(:,kk)   	= WING.FlySTD.GAIN{kk}(:,jj);
        WING.FlyAmpSTD.PhaseDiff  	{jj,1}(:,kk)   	= WING.FlySTD.PhaseDiff{kk}(:,jj);
        
        BODE.FlyAmpSTD.head2wing.GAIN       {jj,1}(:,kk)    = BODE.FlySTD.head2wing.GAIN{kk}(:,jj);
        BODE.FlyAmpSTD.head2wing.PhaseDiff  {jj,1}(:,kk)    = BODE.FlySTD.head2wing.PhaseDiff{kk}(:,jj);
    end
end
clear kk jj
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

    HEAD.GrandMed.Time          = cell2mat((cellfun(@(x) median(x,2),HEAD.FlyAmpMed.Time,'UniformOutput',false))');
    HEAD.GrandMed.Pos           = cell2mat((cellfun(@(x) median(x,2),HEAD.FlyAmpMed.Pos,'UniformOutput',false))');
    HEAD.GrandMed.Vel           = cell2mat((cellfun(@(x) median(x,2),HEAD.FlyAmpMed.Vel,'UniformOutput',false))');
    HEAD.GrandMed.Freq          = cell2mat((cellfun(@(x) median(x,2),HEAD.FlyAmpMed.Freq,'UniformOutput',false))');
    HEAD.GrandMed.Mag           = cell2mat((cellfun(@(x) median(x,2),HEAD.FlyAmpMed.Mag,'UniformOutput',false))');
    HEAD.GrandMed.Phase      	= cell2mat((cellfun(@(x) circ_median(x,2),HEAD.FlyAmpMed.Phase,'UniformOutput',false)))';
    HEAD.GrandMed.COHR.Freq     = cell2mat((cellfun(@(x) median(x,2),HEAD.FlyAmpMed.COHR.Freq,'UniformOutput',false))');
    HEAD.GrandMed.COHR.Mag      = cell2mat((cellfun(@(x) median(x,2),HEAD.FlyAmpMed.COHR.Mag,'UniformOutput',false))');
    HEAD.GrandMed.GAIN      	= cell2mat((cellfun(@(x) median(x,2),HEAD.FlyAmpMed.GAIN,'UniformOutput',false))');
    HEAD.GrandMed.PhaseDiff 	= cell2mat((cellfun(@(x) circ_median(x,2),HEAD.FlyAmpMed.PhaseDiff,'UniformOutput',false)))';
    HEAD.GrandMed.Err.Pos      	= cell2mat((cellfun(@(x) median(x,2),HEAD.FlyAmpMed.Err.Pos,'UniformOutput',false))');
    HEAD.GrandMed.Err.Freq    	= cell2mat((cellfun(@(x) median(x,2),HEAD.FlyAmpMed.Err.Freq ,'UniformOutput',false))');
    HEAD.GrandMed.Err.Mag      	= cell2mat((cellfun(@(x) median(x,2),HEAD.FlyAmpMed.Err.Mag ,'UniformOutput',false))');
    HEAD.GrandMed.Err.Phase    	= cell2mat((cellfun(@(x) circ_median(x,2),HEAD.FlyAmpMed.Err.Phase ,'UniformOutput',false)))';

    WING.GrandMed.Time          = cell2mat((cellfun(@(x) median(x,2),WING.FlyAmpMed.Time,'UniformOutput',false))');
    WING.GrandMed.Pos           = cell2mat((cellfun(@(x) median(x,2),WING.FlyAmpMed.Pos,'UniformOutput',false))');
    WING.GrandMed.Vel           = cell2mat((cellfun(@(x) median(x,2),WING.FlyAmpMed.Vel,'UniformOutput',false))');
    WING.GrandMed.Freq          = cell2mat((cellfun(@(x) median(x,2),WING.FlyAmpMed.Freq,'UniformOutput',false))');
    WING.GrandMed.Mag           = cell2mat((cellfun(@(x) median(x,2),WING.FlyAmpMed.Mag,'UniformOutput',false))');
    WING.GrandMed.Phase      	= cell2mat((cellfun(@(x) circ_median(x,2),WING.FlyAmpMed.Phase,'UniformOutput',false)))';
    WING.GrandMed.COHR.Freq     = cell2mat((cellfun(@(x) median(x,2),WING.FlyAmpMed.COHR.Freq,'UniformOutput',false))');
    WING.GrandMed.COHR.Mag      = cell2mat((cellfun(@(x) median(x,2),WING.FlyAmpMed.COHR.Mag,'UniformOutput',false))');
    WING.GrandMed.GAIN      	= cell2mat((cellfun(@(x) median(x,2),WING.FlyAmpMed.GAIN,'UniformOutput',false))');
    WING.GrandMed.PhaseDiff   	= cell2mat((cellfun(@(x) circ_median(x,2),WING.FlyAmpMed.PhaseDiff,'UniformOutput',false)))';

    BODE.GrandMed.head2wing.GAIN        = cell2mat((cellfun(@(x) median(x,2),BODE.FlyAmpMed.head2wing.GAIN ,'UniformOutput',false))');
    BODE.GrandMed.head2wing.PhaseDiff 	= cell2mat((cellfun(@(x) circ_median(x,2),BODE.FlyAmpMed.head2wing.PhaseDiff ,'UniformOutput',false)))';
% MEAN
%---------------------------------------------------------------------------------------------------------------------------------
    PAT.GrandMean.Time        	= cell2mat((cellfun(@(x) mean(x,2),PAT.FlyAmpMean.Time,'UniformOutput',false))');
    PAT.GrandMean.Pos         	= cell2mat((cellfun(@(x) mean(x,2),PAT.FlyAmpMean.Pos,'UniformOutput',false))');
    PAT.GrandMean.Vel        	= cell2mat((cellfun(@(x) mean(x,2),PAT.FlyAmpMean.Vel,'UniformOutput',false))');
    PAT.GrandMean.Freq       	= cell2mat((cellfun(@(x) mean(x,2),PAT.FlyAmpMean.Freq,'UniformOutput',false))');
    PAT.GrandMean.Mag         	= cell2mat((cellfun(@(x) mean(x,2),PAT.FlyAmpMean.Mag,'UniformOutput',false))');
    PAT.GrandMean.Phase       	= cell2mat((cellfun(@(x) circ_mean(x,[],2),PAT.FlyAmpMean.Phase,'UniformOutput',false)'));

    HEAD.GrandMean.Time       	= cell2mat((cellfun(@(x) mean(x,2),HEAD.FlyAmpMean.Time,'UniformOutput',false))');
    HEAD.GrandMean.Pos          = cell2mat((cellfun(@(x) mean(x,2),HEAD.FlyAmpMean.Pos,'UniformOutput',false))');
    HEAD.GrandMean.Vel      	= cell2mat((cellfun(@(x) mean(x,2),HEAD.FlyAmpMean.Vel,'UniformOutput',false))');
    HEAD.GrandMean.Freq      	= cell2mat((cellfun(@(x) mean(x,2),HEAD.FlyAmpMean.Freq,'UniformOutput',false))');
    HEAD.GrandMean.Mag      	= cell2mat((cellfun(@(x) mean(x,2),HEAD.FlyAmpMean.Mag,'UniformOutput',false))');
    HEAD.GrandMean.Phase     	= cell2mat((cellfun(@(x) circ_mean(x,[],2),HEAD.FlyAmpMean.Phase,'UniformOutput',false)'));
    HEAD.GrandMean.COHR.Freq 	= cell2mat((cellfun(@(x) mean(x,2),HEAD.FlyAmpMean.COHR.Freq,'UniformOutput',false))');
    HEAD.GrandMean.COHR.Mag 	= cell2mat((cellfun(@(x) mean(x,2),HEAD.FlyAmpMean.COHR.Mag,'UniformOutput',false))');
    HEAD.GrandMean.GAIN      	= cell2mat((cellfun(@(x) mean(x,2),HEAD.FlyAmpMean.GAIN,'UniformOutput',false))');
    HEAD.GrandMean.PhaseDiff  	= cell2mat((cellfun(@(x) circ_mean(x,[],2),HEAD.FlyAmpMean.PhaseDiff,'UniformOutput',false)'));
    HEAD.GrandMean.Err.Pos    	= cell2mat((cellfun(@(x) mean(x,2),HEAD.FlyAmpMean.Err.Pos,'UniformOutput',false))');
    HEAD.GrandMean.Err.Freq    	= cell2mat((cellfun(@(x) mean(x,2),HEAD.FlyAmpMean.Err.Freq ,'UniformOutput',false))');
    HEAD.GrandMean.Err.Mag     	= cell2mat((cellfun(@(x) mean(x,2),HEAD.FlyAmpMean.Err.Mag ,'UniformOutput',false))');
    HEAD.GrandMean.Err.Phase   	= cell2mat((cellfun(@(x) circ_mean(x,[],2),HEAD.FlyAmpMean.Err.Phase ,'UniformOutput',false)'));

    WING.GrandMean.Time         = cell2mat((cellfun(@(x) mean(x,2),WING.FlyAmpMean.Time,'UniformOutput',false))');
    WING.GrandMean.Pos      	= cell2mat((cellfun(@(x) mean(x,2),WING.FlyAmpMean.Pos,'UniformOutput',false))');
    WING.GrandMean.Vel          = cell2mat((cellfun(@(x) mean(x,2),WING.FlyAmpMean.Vel,'UniformOutput',false))');
    WING.GrandMean.Freq     	= cell2mat((cellfun(@(x) mean(x,2),WING.FlyAmpMean.Freq,'UniformOutput',false))');
    WING.GrandMean.Mag       	= cell2mat((cellfun(@(x) mean(x,2),WING.FlyAmpMean.Mag,'UniformOutput',false))');
    WING.GrandMean.Phase      	= cell2mat((cellfun(@(x) circ_mean(x,[],2),WING.FlyAmpMean.Phase,'UniformOutput',false)'));
    WING.GrandMean.COHR.Freq  	= cell2mat((cellfun(@(x) mean(x,2),WING.FlyAmpMean.COHR.Freq,'UniformOutput',false))');
    WING.GrandMean.COHR.Mag   	= cell2mat((cellfun(@(x) mean(x,2),WING.FlyAmpMean.COHR.Mag,'UniformOutput',false))');
    WING.GrandMean.GAIN      	= cell2mat((cellfun(@(x) mean(x,2),WING.FlyAmpMean.GAIN,'UniformOutput',false))');
    WING.GrandMean.PhaseDiff  	= cell2mat((cellfun(@(x) circ_mean(x,[],2),WING.FlyAmpMean.PhaseDiff,'UniformOutput',false)'));
  	
    BODE.GrandMean.head2wing.GAIN       = cell2mat((cellfun(@(x) mean(x,2),BODE.FlyAmpMean.head2wing.GAIN ,'UniformOutput',false))');
  	BODE.GrandMean.head2wing.PhaseDiff 	= cell2mat((cellfun(@(x) circ_mean(x,[],2),BODE.FlyAmpMean.head2wing.PhaseDiff ,'UniformOutput',false)')); 
% STD
%---------------------------------------------------------------------------------------------------------------------------------
    PAT.GrandSTD.Time           = cell2mat((cellfun(@(x) std(x,0,2),PAT.FlyAmpSTD.Time,'UniformOutput',false))');
    PAT.GrandSTD.Pos            = cell2mat((cellfun(@(x) std(x,0,2),PAT.FlyAmpSTD.Pos,'UniformOutput',false))');
    PAT.GrandSTD.Vel            = cell2mat((cellfun(@(x) std(x,0,2),PAT.FlyAmpSTD.Vel,'UniformOutput',false))');
    PAT.GrandSTD.Freq           = cell2mat((cellfun(@(x) std(x,0,2),PAT.FlyAmpSTD.Freq,'UniformOutput',false))');
    PAT.GrandSTD.Mag            = cell2mat((cellfun(@(x) std(x,0,2),PAT.FlyAmpSTD.Mag,'UniformOutput',false))');
    PAT.GrandSTD.Phase          = cell2mat((cellfun(@(x) circ_std(x,[],[],2),PAT.FlyAmpSTD.Phase,'UniformOutput',false))');

    HEAD.GrandSTD.Time          = cell2mat((cellfun(@(x) std(x,0,2),HEAD.FlyAmpSTD.Time,'UniformOutput',false))');
    HEAD.GrandSTD.Pos           = cell2mat((cellfun(@(x) std(x,0,2),HEAD.FlyAmpSTD.Pos,'UniformOutput',false))');
    HEAD.GrandSTD.Vel           = cell2mat((cellfun(@(x) std(x,0,2),HEAD.FlyAmpSTD.Vel,'UniformOutput',false))');
    HEAD.GrandSTD.Freq          = cell2mat((cellfun(@(x) std(x,0,2),HEAD.FlyAmpSTD.Freq,'UniformOutput',false))');
    HEAD.GrandSTD.Mag           = cell2mat((cellfun(@(x) std(x,0,2),HEAD.FlyAmpSTD.Mag,'UniformOutput',false))');
    HEAD.GrandSTD.Phase      	= cell2mat((cellfun(@(x) circ_std(x,[],[],2),HEAD.FlyAmpSTD.Phase,'UniformOutput',false))');
    HEAD.GrandSTD.COHR.Freq     = cell2mat((cellfun(@(x) std(x,0,2),HEAD.FlyAmpSTD.COHR.Freq,'UniformOutput',false))');
    HEAD.GrandSTD.COHR.Mag      = cell2mat((cellfun(@(x) std(x,0,2),HEAD.FlyAmpSTD.COHR.Mag,'UniformOutput',false))');
    HEAD.GrandSTD.GAIN      	= cell2mat((cellfun(@(x) std(x,0,2),HEAD.FlyAmpSTD.GAIN,'UniformOutput',false))');
    HEAD.GrandSTD.PhaseDiff  	= cell2mat((cellfun(@(x) circ_std(x,[],[],2),HEAD.FlyAmpSTD.PhaseDiff,'UniformOutput',false))');
    HEAD.GrandSTD.Err.Pos      	= cell2mat((cellfun(@(x) std(x,0,2),HEAD.FlyAmpSTD.Err.Pos,'UniformOutput',false))');
    HEAD.GrandSTD.Err.Freq    	= cell2mat((cellfun(@(x) std(x,0,2),HEAD.FlyAmpSTD.Err.Freq ,'UniformOutput',false))');
    HEAD.GrandSTD.Err.Mag      	= cell2mat((cellfun(@(x) std(x,0,2),HEAD.FlyAmpSTD.Err.Mag ,'UniformOutput',false))');
    HEAD.GrandSTD.Err.Phase    	= cell2mat((cellfun(@(x) circ_std(x,[],[],2),HEAD.FlyAmpSTD.Err.Phase ,'UniformOutput',false))');

    WING.GrandSTD.Time          = cell2mat((cellfun(@(x) std(x,0,2),WING.FlyAmpSTD.Time,'UniformOutput',false))');
    WING.GrandSTD.Pos           = cell2mat((cellfun(@(x) std(x,0,2),WING.FlyAmpSTD.Pos,'UniformOutput',false))');
    WING.GrandSTD.Vel           = cell2mat((cellfun(@(x) std(x,0,2),WING.FlyAmpSTD.Vel,'UniformOutput',false))');
    WING.GrandSTD.Freq          = cell2mat((cellfun(@(x) std(x,0,2),WING.FlyAmpSTD.Freq,'UniformOutput',false))');
    WING.GrandSTD.Mag           = cell2mat((cellfun(@(x) std(x,0,2),WING.FlyAmpSTD.Mag,'UniformOutput',false))');
    WING.GrandSTD.Phase      	= cell2mat((cellfun(@(x) circ_std(x,[],[],2),WING.FlyAmpSTD.Phase,'UniformOutput',false))');
    WING.GrandSTD.COHR.Freq     = cell2mat((cellfun(@(x) std(x,0,2),WING.FlyAmpSTD.COHR.Freq,'UniformOutput',false))');
    WING.GrandSTD.COHR.Mag      = cell2mat((cellfun(@(x) std(x,0,2),WING.FlyAmpSTD.COHR.Mag,'UniformOutput',false))');
    WING.GrandSTD.GAIN      	= cell2mat((cellfun(@(x) std(x,0,2),WING.FlyAmpSTD.GAIN,'UniformOutput',false))');
    WING.GrandSTD.PhaseDiff  	= cell2mat((cellfun(@(x) circ_std(x,[],[],2),WING.FlyAmpSTD.PhaseDiff,'UniformOutput',false))');
    
    BODE.GrandSTD.head2wing.GAIN        = cell2mat((cellfun(@(x) std(x,0,2),BODE.FlyAmpMed.head2wing.GAIN ,'UniformOutput',false))');
    BODE.GrandSTD.head2wing.PhaseDiff 	= cell2mat((cellfun(@(x) circ_std(x,[],[],2),BODE.FlyAmpMed.head2wing.PhaseDiff ,'UniformOutput',false))');
    
%% Save ouputs as structure %%
%---------------------------------------------------------------------------------------------------------------------------------
disp('Saving...')
save([rootdir 'DATA\' filename '.mat'],'PAT','WING','HEAD','BODE','FD','T','n','unq')
disp('DONE')
end