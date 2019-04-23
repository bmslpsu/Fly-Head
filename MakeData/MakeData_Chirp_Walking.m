function [PAT,HEAD,D, N, U] = MakeData_Chirp_Walking(rootdir,filename)
%% MakeData_Chirp_HeadFree: Reads in all raw trials, transforms data, and saves in organized structure for use with figure functions
%   INPUTS:
%       root    : root directory
%   OUTPUTS:
%       PAT     : pattern structure
%       HEAD  	: head structure
%       FD      : file data
%       T       : fly data table
%       n       : field #'s 
%       unq   	: unique fields
%---------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE INPUT %
% rootdir = 'H:\EXPERIMENTS\Experiment_ChirpLog_HeadFree\';
% filename = 'Chirp_Walk_DATA';
%---------------------------------------------------------------------------------------------------------------------------------
%% Setup Directories %%
%---------------------------------------------------------------------------------------------------------------------------------

% Select files
[FILES, PATH] = uigetfile({'*.mat', 'DAQ-files'}, 'Select files', rootdir, 'MultiSelect','on');
FILES = cellstr(FILES)'; % if only one character array >> store in cell

[D,I,N,U] = GetFileData(FILES,'Fly','Trial','Amp');


%% Get Data %%
%---------------------------------------------------------------------------------------------------------------------------------
disp('Loading...')
% Preallocate data cells
HEAD.ALL.Pos = cell(N{1,3},1);
for kk = 1:N{1,1}
    for jj = 1:N{1,3}
        % PATTERN Data
        PAT.Pos{kk,1}{jj,1}         = [];
     	PAT.Time{kk,1}{jj,1}        = [];
        PAT.Vel{kk,1}{jj,1}         = [];
        PAT.VelMean{kk,1}{jj,1} 	= [];
        PAT.VelSTD{kk,1}{jj,1}      = [];
        PAT.Freq{kk,1}{jj,1}        = [];
        PAT.Mag{kk,1}{jj,1}     	= [];
        PAT.Phase{kk,1}{jj,1}       = [];
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
flyFc = 20;
span = 20:1:2100;

for kk = 1:N{1,4}
	filename = fullfile(PATH,FILES{kk}); % full file name
    load(filename,'FlyState','AI','VidTime') % load in fly kinematics & arena voltages
    %-----------------------------------------------------------------------------------------------------------------------------
   
    % Get pattern data %
    pat.Time        = AI{:,1}; % pattern time [s]
    pat.Ts          = mean(diff(pat.Time));
    pat.Fs          = 1/mean(diff(pat.Time)); % pattern sampling frequency [Hz]
    pat.Pos         = 3.75.*round(AI{:,2}.*96/5);
	pat.Vel         = [diff(pat.Pos)/(1/pat.Fs) ; 0]; % pattern velocity [deg/s]
    pat.VelMean     = mean(abs(pat.Vel)); % mean pattern velocity [deg/s]
	pat.VelSTD      = mean(abs(pat.Vel)); % STD pattern velocity [deg/s]
    %-----------------------------------------------------------------------------------------------------------------------------
    % Get head data %
	head.Time       = FlyState{:,1}; % store head time vector [s]
    [b,a]           = butter(2,flyFc/(pat.Fs/2),'low'); % 2nd-order low-pass butterworth filter
    head.Pos        = filtfilt(b,a,FlyState{:,2}); % filter head position
    head.Pos        = head.Pos - mean(head.Pos); % subtract DC component 
    head.Pos        = rad2deg(head.Pos);
    head.Vel        = filtfilt(b,a,[diff(head.Pos)./diff(head.Time) ; 0]); % calculte head vecloity and filter again [deg/s]
    head.VelMean    = mean(abs(head.Vel)); % mean head velocity [deg/s]
    head.VelSTD     = std(abs(head.Vel)); % STD pattern velocity [deg/s]
    %-----------------------------------------------------------------------------------------------------------------------------
	% Interpolate data to match vid.time %
	head.Pos        = interp1(head.Time, head.Pos , VidTime, 'nearest'); % interpolate pattern y-pos to match fly
 	head.Vel        = interp1(head.Time, head.Vel , VidTime, 'nearest'); % interpolate pattern y-pos to match fly
 	pat.Pos         = interp1(pat.Time, pat.Pos  , VidTime, 'nearest'); % interpolate pattern x-pos to match fly
    pat.Vel         = interp1(pat.Time, pat.Vel  , VidTime, 'nearest'); % interpolate pattern x-pos to match fly
   
    head.Time       = VidTime   (span);
    pat.Time        = VidTime   (span);
    head.Pos        = head.Pos  (span);
	head.Vel        = head.Vel  (span);
 	pat.Pos         = pat.Pos   (span);
    pat.Vel         = pat.Vel   (span);
	%-----------------------------------------------------------------------------------------------------------------------------
 	head.Err.Pos    = pat.Pos - head.Pos; % calculate Error between head & pattern (retinal slip) [deg]
    head.Err.Vel    = pat.Vel - head.Vel; % calculate Error between head & pattern (retinal slip) [deg/s]
    %-----------------------------------------------------------------------------------------------------------------------------
    % Convert head, wings, & pattern data to frequency domain using FFT %
    [head.Freq , head.Mag , head.Phase]             = FFT(head.Time,head.Pos);
    [pat.Freq  , pat.Mag  , pat.Phase ]             = FFT(pat.Time,pat.Pos);
	[head.Err.Freq , head.Err.Mag, head.Err.Phase]  = FFT(head.Time,head.Err.Pos);
    %-----------------------------------------------------------------------------------------------------------------------------
    % Calculate coherence %
    [head.cohr.mag,head.cohr.f] = mscohere(pat.Pos , head.Pos ,[],[] , head.Freq , pat.Fs);
	
    %-----------------------------------------------------------------------------------------------------------------------------
	% Calculate BODE gain & phase difference for head & wings %
    head.GAIN   = medfilt1(head.Mag./pat.Mag,5);
    head.PHASE  = medfilt1(-(pat.Phase - head.Phase),5);
   
    %-----------------------------------------------------------------------------------------------------------------------------
    % Store data in cells %
  	% PATTERN
	PAT.Time        {I{kk,1},1}{I{kk,3},1}(:,end+1) = pat.Time;
	PAT.Pos         {I{kk,1},1}{I{kk,3},1}(:,end+1) = pat.Pos;
	PAT.Vel         {I{kk,1},1}{I{kk,3},1}(:,end+1) = pat.Vel;
   	PAT.VelMean     {I{kk,1},1}{I{kk,3},1}(:,end+1) = pat.VelMean;
	PAT.VelSTD      {I{kk,1},1}{I{kk,3},1}(:,end+1) = pat.VelSTD;
    PAT.Freq        {I{kk,1},1}{I{kk,3},1}(:,end+1) = pat.Freq;
	PAT.Mag         {I{kk,1},1}{I{kk,3},1}(:,end+1) = pat.Mag;
	PAT.Phase       {I{kk,1},1}{I{kk,3},1}(:,end+1) = pat.Phase;
    
    % HEAD
    HEAD.ALL.Pos    {I{kk,1}}(:,end+1) = head.Pos;
	HEAD.Time       {I{kk,1},1}{I{kk,3},1}(:,end+1) = head.Time;
	HEAD.Pos        {I{kk,1},1}{I{kk,3},1}(:,end+1) = head.Pos;
    HEAD.Vel        {I{kk,1},1}{I{kk,3},1}(:,end+1) = head.Vel;
    HEAD.VelMean    {I{kk,1},1}{I{kk,3},1}(:,end+1) = head.VelMean;
    HEAD.VelSTD     {I{kk,1},1}{I{kk,3},1}(:,end+1) = head.VelSTD;
	HEAD.Err.Pos    {I{kk,1},1}{I{kk,3},1}(:,end+1) = head.Err.Pos;
    HEAD.Err.Freq 	{I{kk,1},1}{I{kk,3},1}(:,end+1) = head.Err.Freq;
    HEAD.Err.Mag    {I{kk,1},1}{I{kk,3},1}(:,end+1) = head.Err.Mag;
 	HEAD.Err.Phase 	{I{kk,1},1}{I{kk,3},1}(:,end+1) = head.Err.Phase;
	HEAD.Freq       {I{kk,1},1}{I{kk,3},1}(:,end+1) = head.Freq;
	HEAD.Mag        {I{kk,1},1}{I{kk,3},1}(:,end+1) = head.Mag;
	HEAD.Phase      {I{kk,1},1}{I{kk,3},1}(:,end+1) = head.Phase;
	HEAD.COHR.Freq  {I{kk,1},1}{I{kk,3},1}(:,end+1) = head.cohr.f;
    HEAD.COHR.Mag   {I{kk,1},1}{I{kk,3},1}(:,end+1) = head.cohr.mag;
    HEAD.GAIN       {I{kk,1},1}{I{kk,3},1}(:,end+1) = head.GAIN;
	HEAD.PHASE      {I{kk,1},1}{I{kk,3},1}(:,end+1) = head.PHASE;
   
end
clear jj kk a b VidTime FlyState AI head pat 
disp('LOADING DONE')
%% FLY Stats by Fly %%
%---------------------------------------------------------------------------------------------------------------------------------
for kk = 1:N{1,1}
    for jj = 1:N{1,3}
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
        HEAD.FlyMed.GAIN        {kk,1}(:,jj)	= median(HEAD.GAIN{kk}{jj},2);
        HEAD.FlyMed.PHASE      	{kk,1}(:,jj)	= circ_median(HEAD.PHASE{kk}{jj},2)';
        HEAD.FlyMed.Err.Pos    	{kk,1}(:,jj)	= median(HEAD.Err.Pos{kk}{jj},2);
        HEAD.FlyMed.Err.Freq  	{kk,1}(:,jj)	= median(HEAD.Err.Freq{kk}{jj},2);
        HEAD.FlyMed.Err.Mag    	{kk,1}(:,jj)	= median(HEAD.Err.Mag{kk}{jj},2);
    	HEAD.FlyMed.Err.Phase 	{kk,1}(:,jj)	= circ_median(HEAD.Err.Phase{kk}{jj},2)';
   
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
        HEAD.FlyMean.PHASE      {kk,1}(:,jj)	= circ_mean(HEAD.PHASE{kk}{jj},[],2)';
        HEAD.FlyMean.Err.Pos 	{kk,1}(:,jj)	= mean(HEAD.Err.Pos{kk}{jj},2);
        HEAD.FlyMean.Err.Freq  	{kk,1}(:,jj)	= mean(HEAD.Err.Freq{kk}{jj},2);
        HEAD.FlyMean.Err.Mag   	{kk,1}(:,jj)	= mean(HEAD.Err.Mag{kk}{jj},2);
    	HEAD.FlyMean.Err.Phase 	{kk,1}(:,jj)	= circ_mean(HEAD.Err.Phase{kk}{jj},[],2)'; 
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
        HEAD.FlySTD.PHASE      	{kk,1}(:,jj)	= circ_std(HEAD.PHASE{kk}{jj},[],[],2);

        HEAD.FlySTD.Err.Pos    	{kk,1}(:,jj)	= std(HEAD.Err.Pos{kk}{jj},0,2);
        HEAD.FlySTD.Err.Freq  	{kk,1}(:,jj)	= std(HEAD.Err.Freq{kk}{jj},0,2);
        HEAD.FlySTD.Err.Mag    	{kk,1}(:,jj)	= std(HEAD.Err.Mag{kk}{jj},0,2);
    	HEAD.FlySTD.Err.Phase 	{kk,1}(:,jj)	= circ_std(HEAD.Err.Phase{kk}{jj},[],[],2);          
    end
end
clear kk jj
%% FLY stats by Amplitude %%
%---------------------------------------------------------------------------------------------------------------------------------
for jj = 1:N{1,3}
    for kk = 1:N{1,1}
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
        HEAD.FlyAmpMed.PHASE        {jj,1}(:,kk)   	= HEAD.FlyMed.PHASE{kk}(:,jj);
        HEAD.FlyAmpMed.Err.Pos      {jj,1}(:,kk) 	= HEAD.FlyMed.Err.Pos{kk}(:,jj);
    	HEAD.FlyAmpMed.Err.Freq     {jj,1}(:,kk) 	= HEAD.FlyMed.Err.Freq{kk}(:,jj);
    	HEAD.FlyAmpMed.Err.Mag      {jj,1}(:,kk) 	= HEAD.FlyMed.Err.Mag{kk}(:,jj);
    	HEAD.FlyAmpMed.Err.Phase    {jj,1}(:,kk) 	= HEAD.FlyMed.Err.Phase{kk}(:,jj);
        
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
        HEAD.FlyAmpMean.PHASE     	{jj,1}(:,kk)   	= HEAD.FlyMean.PHASE{kk}(:,jj);
        HEAD.FlyAmpMean.Err.Pos   	{jj,1}(:,kk) 	= HEAD.FlyMean.Err.Pos{kk}(:,jj);
    	HEAD.FlyAmpMean.Err.Freq  	{jj,1}(:,kk) 	= HEAD.FlyMean.Err.Freq{kk}(:,jj);
    	HEAD.FlyAmpMean.Err.Mag 	{jj,1}(:,kk) 	= HEAD.FlyMean.Err.Mag{kk}(:,jj);
    	HEAD.FlyAmpMean.Err.Phase  	{jj,1}(:,kk) 	= HEAD.FlyMean.Err.Phase{kk}(:,jj);

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
        HEAD.FlyAmpSTD.PHASE     	{jj,1}(:,kk)   	= HEAD.FlySTD.PHASE{kk}(:,jj);
        HEAD.FlyAmpSTD.Err.Pos   	{jj,1}(:,kk) 	= HEAD.FlySTD.Err.Pos{kk}(:,jj);
    	HEAD.FlyAmpSTD.Err.Freq  	{jj,1}(:,kk) 	= HEAD.FlySTD.Err.Freq{kk}(:,jj);
    	HEAD.FlyAmpSTD.Err.Mag      {jj,1}(:,kk) 	= HEAD.FlySTD.Err.Mag{kk}(:,jj);
    	HEAD.FlyAmpSTD.Err.Phase  	{jj,1}(:,kk) 	= HEAD.FlySTD.Err.Phase{kk}(:,jj);
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
    PAT.GrandMed.Phase          = cell2mat((cellfun(@(x) median(x,2),PAT.FlyAmpMed.Phase,'UniformOutput',false)))';

    HEAD.GrandMed.Time          = cell2mat((cellfun(@(x) median(x,2),HEAD.FlyAmpMed.Time,'UniformOutput',false))');
    HEAD.GrandMed.Pos           = cell2mat((cellfun(@(x) median(x,2),HEAD.FlyAmpMed.Pos,'UniformOutput',false))');
    HEAD.GrandMed.Vel           = cell2mat((cellfun(@(x) median(x,2),HEAD.FlyAmpMed.Vel,'UniformOutput',false))');
    HEAD.GrandMed.Freq          = cell2mat((cellfun(@(x) median(x,2),HEAD.FlyAmpMed.Freq,'UniformOutput',false))');
    HEAD.GrandMed.Mag           = cell2mat((cellfun(@(x) median(x,2),HEAD.FlyAmpMed.Mag,'UniformOutput',false))');
    HEAD.GrandMed.Phase      	= cell2mat((cellfun(@(x) median(x,2),HEAD.FlyAmpMed.Phase,'UniformOutput',false)))';
    HEAD.GrandMed.COHR.Freq     = cell2mat((cellfun(@(x) median(x,2),HEAD.FlyAmpMed.COHR.Freq,'UniformOutput',false))');
    HEAD.GrandMed.COHR.Mag      = cell2mat((cellfun(@(x) median(x,2),HEAD.FlyAmpMed.COHR.Mag,'UniformOutput',false))');
    HEAD.GrandMed.GAIN      	= cell2mat((cellfun(@(x) median(x,2),HEAD.FlyAmpMed.GAIN,'UniformOutput',false))');
    HEAD.GrandMed.PHASE      	= cell2mat((cellfun(@(x) median(x,2),HEAD.FlyAmpMed.PHASE,'UniformOutput',false)))';
    HEAD.GrandMed.Err.Pos      	= cell2mat((cellfun(@(x) median(x,2),HEAD.FlyAmpMed.PHASE,'UniformOutput',false))');
    HEAD.GrandMed.Err.Freq    	= cell2mat((cellfun(@(x) median(x,2),HEAD.FlyAmpMed.Err.Freq ,'UniformOutput',false))');
    HEAD.GrandMed.Err.Mag      	= cell2mat((cellfun(@(x) median(x,2),HEAD.FlyAmpMed.Err.Mag ,'UniformOutput',false))');
    HEAD.GrandMed.Err.Phase    	= cell2mat((cellfun(@(x) median(x,2),HEAD.FlyAmpMed.Err.Phase ,'UniformOutput',false)))';

% MEAN
%---------------------------------------------------------------------------------------------------------------------------------
    PAT.GrandMean.Time        	= cell2mat((cellfun(@(x) mean(x,2),PAT.FlyAmpMean.Time,'UniformOutput',false))');
    PAT.GrandMean.Pos         	= cell2mat((cellfun(@(x) mean(x,2),PAT.FlyAmpMean.Pos,'UniformOutput',false))');
    PAT.GrandMean.Vel        	= cell2mat((cellfun(@(x) mean(x,2),PAT.FlyAmpMean.Vel,'UniformOutput',false))');
    PAT.GrandMean.Freq       	= cell2mat((cellfun(@(x) mean(x,2),PAT.FlyAmpMean.Freq,'UniformOutput',false))');
    PAT.GrandMean.Mag         	= cell2mat((cellfun(@(x) mean(x,2),PAT.FlyAmpMean.Mag,'UniformOutput',false))');
    PAT.GrandMean.Phase       	= cell2mat((cellfun(@(x) mean(x,2),PAT.FlyAmpMean.Phase,'UniformOutput',false)'));

    HEAD.GrandMean.Time       	= cell2mat((cellfun(@(x) mean(x,2),HEAD.FlyAmpMean.Time,'UniformOutput',false))');
    HEAD.GrandMean.Pos          = cell2mat((cellfun(@(x) mean(x,2),HEAD.FlyAmpMean.Pos,'UniformOutput',false))');
    HEAD.GrandMean.Vel      	= cell2mat((cellfun(@(x) mean(x,2),HEAD.FlyAmpMean.Vel,'UniformOutput',false))');
    HEAD.GrandMean.Freq      	= cell2mat((cellfun(@(x) mean(x,2),HEAD.FlyAmpMean.Freq,'UniformOutput',false))');
    HEAD.GrandMean.Mag      	= cell2mat((cellfun(@(x) mean(x,2),HEAD.FlyAmpMean.Mag,'UniformOutput',false))');
    HEAD.GrandMean.Phase     	= cell2mat((cellfun(@(x) mean(x,2),HEAD.FlyAmpMean.Phase,'UniformOutput',false)'));
    HEAD.GrandMean.COHR.Freq 	= cell2mat((cellfun(@(x) mean(x,2),HEAD.FlyAmpMean.COHR.Freq,'UniformOutput',false))');
    HEAD.GrandMean.COHR.Mag 	= cell2mat((cellfun(@(x) mean(x,2),HEAD.FlyAmpMean.COHR.Mag,'UniformOutput',false))');
    HEAD.GrandMean.GAIN      	= cell2mat((cellfun(@(x) mean(x,2),HEAD.FlyAmpMean.GAIN,'UniformOutput',false))');
    HEAD.GrandMean.PHASE      	= cell2mat((cellfun(@(x) mean(x,2),HEAD.FlyAmpMean.PHASE,'UniformOutput',false)'));
    HEAD.GrandMean.Err.Pos    	= cell2mat((cellfun(@(x) mean(x,2),HEAD.FlyAmpMean.PHASE,'UniformOutput',false))');
    HEAD.GrandMean.Err.Freq    	= cell2mat((cellfun(@(x) mean(x,2),HEAD.FlyAmpMean.Err.Freq ,'UniformOutput',false))');
    HEAD.GrandMean.Err.Mag     	= cell2mat((cellfun(@(x) mean(x,2),HEAD.FlyAmpMean.Err.Mag ,'UniformOutput',false))');
    HEAD.GrandMean.Err.Phase   	= cell2mat((cellfun(@(x) mean(x,2),HEAD.FlyAmpMean.Err.Phase ,'UniformOutput',false)'));

   
% STD
%---------------------------------------------------------------------------------------------------------------------------------
    PAT.GrandSTD.Time           = cell2mat((cellfun(@(x) std(x,0,2),PAT.FlyAmpMed.Time,'UniformOutput',false))');
    PAT.GrandSTD.Pos            = cell2mat((cellfun(@(x) std(x,0,2),PAT.FlyAmpMed.Pos,'UniformOutput',false))');
    PAT.GrandSTD.Vel            = cell2mat((cellfun(@(x) std(x,0,2),PAT.FlyAmpMed.Vel,'UniformOutput',false))');
    PAT.GrandSTD.Freq           = cell2mat((cellfun(@(x) std(x,0,2),PAT.FlyAmpMed.Freq,'UniformOutput',false))');
    PAT.GrandSTD.Mag            = cell2mat((cellfun(@(x) std(x,0,2),PAT.FlyAmpMed.Mag,'UniformOutput',false))');
    PAT.GrandSTD.Phase          = cell2mat((cellfun(@(x) std(x,0,2),PAT.FlyAmpMed.Phase,'UniformOutput',false))');

    HEAD.GrandSTD.Time          = cell2mat((cellfun(@(x) std(x,0,2),HEAD.FlyAmpMed.Time,'UniformOutput',false))');
    HEAD.GrandSTD.Pos           = cell2mat((cellfun(@(x) std(x,0,2),HEAD.FlyAmpMed.Pos,'UniformOutput',false))');
    HEAD.GrandSTD.Vel           = cell2mat((cellfun(@(x) std(x,0,2),HEAD.FlyAmpMed.Vel,'UniformOutput',false))');
    HEAD.GrandSTD.Freq          = cell2mat((cellfun(@(x) std(x,0,2),HEAD.FlyAmpMed.Freq,'UniformOutput',false))');
    HEAD.GrandSTD.Mag           = cell2mat((cellfun(@(x) std(x,0,2),HEAD.FlyAmpMed.Mag,'UniformOutput',false))');
    HEAD.GrandSTD.Phase      	= cell2mat((cellfun(@(x) std(x,0,2),HEAD.FlyAmpMed.Phase,'UniformOutput',false))');
    HEAD.GrandSTD.COHR.Freq     = cell2mat((cellfun(@(x) std(x,0,2),HEAD.FlyAmpMed.COHR.Freq,'UniformOutput',false))');
    HEAD.GrandSTD.COHR.Mag      = cell2mat((cellfun(@(x) std(x,0,2),HEAD.FlyAmpMed.COHR.Mag,'UniformOutput',false))');
    HEAD.GrandSTD.GAIN      	= cell2mat((cellfun(@(x) std(x,0,2),HEAD.FlyAmpMed.GAIN,'UniformOutput',false))');
    HEAD.GrandSTD.PHASE      	= cell2mat((cellfun(@(x) std(x,0,2),HEAD.FlyAmpMed.PHASE,'UniformOutput',false))');
    HEAD.GrandSTD.Err.Pos      	= cell2mat((cellfun(@(x) std(x,0,2),HEAD.FlyAmpMed.PHASE,'UniformOutput',false))');
    HEAD.GrandSTD.Err.Freq    	= cell2mat((cellfun(@(x) std(x,0,2),HEAD.FlyAmpMed.Err.Freq ,'UniformOutput',false))');
    HEAD.GrandSTD.Err.Mag      	= cell2mat((cellfun(@(x) std(x,0,2),HEAD.FlyAmpMed.Err.Mag ,'UniformOutput',false))');
    HEAD.GrandSTD.Err.Phase    	= cell2mat((cellfun(@(x) std(x,0,2),HEAD.FlyAmpMed.Err.Phase ,'UniformOutput',false))');

%% Save ouputs as structure %%
%---------------------------------------------------------------------------------------------------------------------------------
disp('Saving...')
save([rootdir 'DATA\' filename '.mat'],'PAT','HEAD','D','N','U')
disp('SAVING DONE')

end