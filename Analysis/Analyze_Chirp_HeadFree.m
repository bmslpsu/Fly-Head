function [] = Analyze_Chirp_HeadFree()
%% Analyze Chirp Data:
%   INPUTS:
%       root: root directory
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
showplot.Time  = 1;
showplot.Freq  = 1;
showplot.Bode  = 0;
showplot.Coher = 0;
showplot.FFT   = 0;
%% Setup Directories %%
%---------------------------------------------------------------------------------------------------------------------------------
root.pat = 'E:\Experiment_HeadExcitation\Chirp\HeadFree\';
root.head = [root.pat '\Vid\Angles\'];

% Select files
[FILES, PATH.head] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root.head, 'MultiSelect','on');
FILES = FILES';

PATH.pat = uigetdir(root.pat);
PATH.pat = [PATH.pat '\'];

%% Process File Data %%
%---------------------------------------------------------------------------------------------------------------------------------
% Read in data from head file names: [Fly=FLy#, Trial=trial# for each fly]nTrial = length(FILES);     % total # of trials
nTrial = length(FILES);     % total # of trials
Fly = zeros(nTrial,1); Trial = zeros(nTrial,1); Amp = zeros(nTrial,1); HEAD.FileCells = cell(nTrial,6);% preallocate arrays
for jj = 1:nTrial
    temp = textscan(char(FILES{jj}), '%s', 'delimiter', '_.'); temp = temp{1} ; % read individual strings into temp variable
    HEAD.FileCells(jj,:) = {temp{1} temp{2} temp{3} temp{4} temp{5} temp{6}}; % separate strings into six rows with one item in each cell
    Fly(jj,1) = str2double(temp{2}); % store fly #
    Trial(jj,1) = str2double(temp{4}); % store trial #
    
    if ~isnan(str2double(temp{7}))
        Amp(jj,1) = str2double([temp{6} '.' temp{7}]); % store frequency
    else
        Amp(jj,1) = str2double(temp{6});
    end
end
%% Set up indexing convention for files %%
% Normalize to start at fly#1 and increment by 1 for each fy (same for trials)
%---------------------------------------------------------------------------------------------------------------------------------
uFly = sort(unique(Fly)); % original # of each unique fly
nFly = length(uFly); % # flies
trialFly = cell(nFly,2);
for kk = 1:nFly
    trialFly{kk,1} = uFly(kk);  % fly # label
	trialFly{kk,2} = length(find(Fly==uFly(kk))); % # trials per fly
end
% Make indexing array for flies
newFly = 1:nFly;
idxFly = zeros(length(Fly),1);
pp = 1;
for kk = 1:nFly
    idxFly(pp:pp+trialFly{kk,2}-1,1) = newFly(kk); % fly index
    pp = pp+trialFly{kk,2};
end
% Make indexing array for trials
pp = 1;
idxTrial = zeros(nTrial,1);
for kk = 1:nFly
   idxTrial(pp:pp+trialFly{kk,2}-1) = 1:trialFly{kk,2}; % trial index
   pp = pp+trialFly{kk,2};
end
% Make indexing array for frequencies
uAmp = sort(unique(Amp)); % find all unique frequencies
nAmp = length(uAmp); % # of unique frequencies
idxAmp = zeros(nTrial,1);
for kk = 1:nAmp
   idx = find(Amp == uAmp(kk));
   idxAmp(idx) = kk; % frequency index
end

fprintf('Total Flies''s: %i \n',nFly) ; fprintf('Total Trials''s: %i \n',nTrial)
T = cell2table(trialFly,'VariableNames',{'Fly','Trials'});
disp(T)
%% Get Data %%
%---------------------------------------------------------------------------------------------------------------------------------
disp('Loading Data...')
clear HEAD PAT WING head pat wings Amplitude
SI = 1;
EI = 241;
Amplitude = zeros(nTrial,1);
% Preallocate data cells %
for kk = 1:nFly % # of flys
    for jj = 1:nAmp % # of amplitudes
        HEAD.Time{kk,1}{jj,1}       = [];
        HEAD.Pos{kk,1}{jj,1}        = [];
        HEAD.Vel{kk,1}{jj,1}        = [];
        HEAD.VelMean{kk,1}{jj,1}    = [];
        HEAD.VelSTD{kk,1}{jj,1}     = [];
        HEAD.IntErr{kk,1}{jj,1}     = [];
        HEAD.Freq{kk,1}{jj,1}       = [];
        HEAD.Mag{kk,1}{jj,1}      	= [];
        HEAD.Phase{kk,1}{jj,1}      = [];
        HEAD.GAIN{kk,1}{jj,1}       = [];
        HEAD.PHASE{kk,1}{jj,1}      = [];
        HEAD.Err.Pos{kk,1}{jj,1}	= [];
        HEAD.Err.Freq{kk,1}{jj,1}	= [];
        HEAD.Err.Mag{kk,1}{jj,1}	= [];
        HEAD.Err.Phase{kk,1}{jj,1}	= [];
        HEAD.COHR.Freq{kk,1}{jj,1}  = [];
        HEAD.COHR.Mag{kk,1}{jj,1}   = [];
        
        PAT.Time{kk,1}{jj,1}        = [];
        PAT.Pos{kk,1}{jj,1}         = [];
        PAT.Freq{kk,1}{jj,1}        = [];
        PAT.Mag{kk,1}{jj,1}     	= [];
        PAT.Phase{kk,1}{jj,1}       = [];
        PAT.Vel{kk,1}{jj,1}         = [];
        
        WING.Time{kk,1}{jj,1}       = [];
        WING.Pos{kk,1}{jj,1}        = [];
        WING.Freq{kk,1}{jj,1}       = [];
        WING.Mag{kk,1}{jj,1}     	= [];
        WING.Phase{kk,1}{jj,1}      = [];
        WING.GAIN{kk,1}{jj,1}       = [];
        WING.PHASE{kk,1}{jj,1}      = [];
        WING.Err.GAIN{kk,1}{jj,1} 	= [];
        WING.Err.PHASE{kk,1}{jj,1} 	= [];
        WING.COHR.Freq{kk,1}{jj,1}  = [];
        WING.COHR.Mag{kk,1}{jj,1}   = [];
    end
end
% Save time & frequency domain data in cells %
for kk = 1:nTrial
    clear head pat func data t_v t_p sineData hAngles % clear temporary variables
    %-----------------------------------------------------------------------------------------------------------------------------
    % Load head & DAQ data %
    load([PATH.head  FILES{kk}]); % load head angles % time arrays
	load([PATH.pat   FILES{kk}]); % load pattern x-position
	%-----------------------------------------------------------------------------------------------------------------------------
	% Get head data %
	head.Time = t_v; % store head time vector
    head.Fs = 1/mean(diff(head.Time)); % sampling frequency [Hz]
    head.Fc = 20; % cutoff frequency [Hz]
    [b,a] = butter(2,head.Fc/(head.Fs/2)); % butterworth filter
    head.Pos = filtfilt(b,a,hAngles); % filter head position
    head.Pos = head.Pos - mean(head.Pos); % subtract DC component
    head.Vel = filtfilt(b,a,[diff(head.Pos)./diff(head.Time) ; 0]);
    head.VelMean = mean(abs(head.Vel));
    head.VelSTD = std(abs(head.Vel));
  	%-----------------------------------------------------------------------------------------------------------------------------
    % Get wing data from DAQ %
  	K = 1; % WBA factor
    wing.Full.Time  = t_p; % wing time 
    wing.Fs = 1/mean(diff(wing.Full.Time)); % sampling frequency [Hz]
    wing.Fc = 20; % cutoff frequency [Hz]
    [b,a] = butter(2,wing.Fc/(wing.Fs/2)); % butterworth filter
	wing.Left       = filtfilt(b,a,(data(:,4))); % left wing
    wing.Right      = filtfilt(b,a,(data(:,5))); % right wing
    wing.Full.Pos 	= (K*(wing.Left - wing.Right)); % WBA (L-R)
  	wing.Full.Pos   = wing.Full.Pos - mean(wing.Full.Pos); % subtract mean
  	wing.Full.Vel   = [diff(wing.Full.Pos)./diff(wing.Full.Time) ; 0];
    
    interval = floor(length(wing.Full.Time)/length(head.Time));
    wing.Pos  = wing.Full.Pos(1:interval:end-1);
    wing.Vel  = wing.Full.Vel(1:interval:end-1);
    wing.Time = wing.Full.Time(1:interval:end-1);
           
 	wing.f = 100*(data(:,6)); % wing beat frequency
    if min(wing.f)<180
        fprintf('Low WBF: Fly %i Trial %i \n',Fly(kk),Trial(kk))
    end
 	%-----------------------------------------------------------------------------------------------------------------------------
    % Get pattern data from DAQ %
    pat.Full.Time = t_p;
    pat.Full.Pos  = 3.75*round((96/10)*(data(:,2)' - mean(data(:,2))))'; % get pattern x-pos: subtract mean and convert to deg
	pat.Full.Vel = [diff(pat.Full.Pos)./diff(pat.Full.Time) ; 0];

    interval = floor(length(pat.Full.Time)/length(head.Time));
    pat.Pos  = pat.Full.Pos(1:interval:end-1);
	pat.Time = pat.Full.Time(1:interval:end-1);
	pat.Vel = [diff(pat.Pos)./diff(pat.Time) ; 0];
	Amplitude(kk,1) = max(pat.Pos);
    %-----------------------------------------------------------------------------------------------------------------------------
    % Convert head, wings, & pattern data to frequency domain %
    [head.Freq , head.Mag, head.Phase]   = FFT(head.Time,head.Pos);
    [wing.Freq , wing.Mag , wing.Phase]  = FFT(wing.Time,wing.Pos);
    [pat.Freq  , pat.Mag  , pat.Phase ]  = FFT(pat.Time,pat.Pos);
    
    head.Freq   = head.Freq(SI:EI);
    head.Mag    = head.Mag(SI:EI);
    head.Phase  = head.Phase(SI:EI);
    
    wing.Freq   = wing.Freq(SI:EI);
    wing.Mag    = wing.Mag(SI:EI);
    wing.Phase  = wing.Phase(SI:EI);
	
    pat.Freq   = pat.Freq(SI:EI);
    pat.Mag    = pat.Mag(SI:EI);
    pat.Phase  = pat.Phase(SI:EI);
	%-----------------------------------------------------------------------------------------------------------------------------
    % Calculate BODE gain % phase %   
    head.bode.mag = medfilt1(head.Mag./pat.Mag,1); % head gain
    head.bode.phase = -1*medfilt1(pat.Phase - head.Phase,1); % head phase
    %-----------------------------------------------------------------------------------------------------------------------------
    % Calculate coherence %
    [head.cohr.mag,head.cohr.f] = mscohere(pat.Pos ,head.Pos,[],[],head.Freq,head.Fs);
	[wing.cohr.mag,wing.cohr.f] = mscohere(pat.Pos ,wing.Pos ,[],[],wing.Freq ,head.Fs);
    %-----------------------------------------------------------------------------------------------------------------------------
    % Calculate Error to Wing BODE after removing head movements %
    head.Err.Pos = pat.Pos - head.Pos;
    [head.Err.Freq , head.Err.Mag, head.Err.Phase] = FFT(head.Time,head.Err.Pos);
    
    head.Err.Freq  = head.Err.Freq(SI:EI);
   	head.Err.Mag   = head.Err.Mag(SI:EI);
    head.Err.Phase = head.Err.Phase(SI:EI);

    wing.Err.bode.mag   = medfilt1(wing.Mag./head.Err.Mag,30);
    wing.Err.bode.phase = -1*medfilt1(head.Err.Phase - wing.Phase,30);
    %-----------------------------------------------------------------------------------------------------------------------------
    % Store data in cells %
    % Head
	HEAD.Time     {idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = head.Time;
	HEAD.Pos  	  {idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = head.Pos;
    HEAD.Vel  	  {idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = head.Vel;
    HEAD.VelMean  {idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = head.VelMean;
    HEAD.VelSTD   {idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = head.VelSTD;

	HEAD.Freq   {idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = head.Freq;
	HEAD.Mag 	{idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = head.Mag;
	HEAD.Phase	{idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = head.Phase;
    
    HEAD.GAIN   {idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = head.bode.mag;
	HEAD.PHASE  {idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = head.bode.phase;
        
    HEAD.Err.Pos    {idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = head.Err.Pos;
    HEAD.Err.Freq 	{idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = head.Err.Freq;
    HEAD.Err.Mag    {idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = head.Err.Mag;
 	HEAD.Err.Phase 	{idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = head.Err.Phase;

    HEAD.COHR.Freq  {idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = head.cohr.f;
    HEAD.COHR.Mag   {idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = head.cohr.mag;
    
    % Pattern
	PAT.Time  	{idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = pat.Full.Time;
	PAT.Pos    	{idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = pat.Pos;
	PAT.Vel    	{idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = pat.Vel;

    PAT.Freq  	{idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = pat.Freq;
	PAT.Mag    	{idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = pat.Mag;
	PAT.Phase	{idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = pat.Phase;
    
    % Wings
	WING.Pos	{idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = wing.Pos;
	WING.Time	{idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = wing.Time;

	WING.Freq  {idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = wing.Freq;
	WING.Mag    {idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = wing.Mag;
	WING.Phase {idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = wing.Phase;
    
    WING.Err.GAIN   {idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = wing.Err.bode.mag;
  	WING.Err.PHASE  {idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = wing.Err.bode.phase;

	WING.COHR.Freq  {idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = wing.cohr.f;
    WING.COHR.Mag   {idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = wing.cohr.mag;    
    %-----------------------------------------------------------------------------------------------------------------------------
    colmn = 4;
    if showplot.Time
        figure (100)
            subplot(ceil(nTrial/colmn),colmn,kk) ; hold on
                title(['Fly ' num2str(Fly(kk)) ' Amp ' num2str(Amp(kk))])
                plot(pat.Time,pat.Pos,'k')
                plot(head.Time,head.Pos,'b')
    end
                plot(wing.Time,wing.Pos,'r')
                box on
    if showplot.Freq
        figure (104)
            subplot(ceil(nTrial/colmn),colmn,kk) ; hold on
                title(['Fly ' num2str(Fly(kk)) ' Amp ' num2str(Amp(kk))])
                plot(pat.Freq,  pat.Mag ,'k')
                plot(head.Freq, head.Mag,'b')
             	plot(wing.Freq, wing.Mag,'r')
                xlim([0.5 11.5])
                box on
        figure (105)
            subplot(ceil(nTrial/colmn),colmn,kk) ; hold on
                title(['Fly ' num2str(Fly(kk)) ' Amp ' num2str(Amp(kk))])
                plot(pat.Freq,  pat.Phase ,'k')
                plot(head.Freq, head.Phase ,'b')
                plot(wing.Freq, wing.Phase,'r')
                xlim([0.5 11.5])
                box on    
    end
    
    if showplot.Bode
        figure (109)
            subplot(ceil(nTrial/colmn),colmn,kk) ; hold on
                title(['Fly ' num2str(Fly(kk)) ' Amp ' num2str(Amp(kk))])
                plot(head.Freq, head.bode.mag ,'k')
                xlim([0.5 11.5])
                box on
        figure (110)
            subplot(ceil(nTrial/colmn),colmn,kk) ; hold on
                title(['Fly ' num2str(Fly(kk)) ' Amp ' num2str(Amp(kk))])
                plot(head.Freq, head.bode.phase ,'k')
                xlim([0.5 11.5])
                box on
    end
    %-----------------------------------------------------------------------------------------------------------------------------
end
disp('DONE')

%% Example Trial %%
%---------------------------------------------------------------------------------------------------------------------------------
ttt = 0:(1/500):20;
Stim = 15*chirp(ttt,0.1,20,12,'Logarithmic',45);
ff = 4;
HH = HEAD.Pos{ff}{3}(:,2);
HT = HEAD.Time{ff}{3}(:,2);
WT = WING.Pos{ff}{3}(:,2);

figure (322) ; clf ; hold on ; box on
set(gcf,'Position',[100 100 900 600])
xlabel('Time (s)')
xticks([0 20])
set(gca,'fontsize',28)

yyaxis left ; set(gca,'YColor','b') ; ylabel(['Angle (' char(176) ')']) ; ylim([-15 15]) ; yticks([-15 0 15])
    plot(ttt,0*Stim,'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',2)
    hh.Stim = plot(ttt,Stim,'-k','LineWidth',2,'HandleVisibility','off');
    hh.Head =  plot(HT,HH,'-b','LineWidth',3);
    
yyaxis right ; set(gca,'YColor','r') ; ylabel('Voltage (V)') ; ylim([-6 6]) ; yticks([-6 0 6])
    hh.Wing = plot(HT,WT,'Color','r','LineWidth',3);

legend([hh.Stim,hh.Head,hh.Wing],{['Visual Scene (' char(176) ')'],['Head (' char(176) ')'],'Wings (V)'},...
    'Orientation','horizontal','Box','off','Location','north')

%% Averages %%
%---------------------------------------------------------------------------------------------------------------------------------
clc
for kk = 1:nFly
    for jj = 1:nAmp
     	PAT.FlyMed.Freq         {kk,1}   	= zeros(length(PAT.Freq{kk}{jj}),nAmp);
        PAT.FlyMed.Phase        {kk,1}   	= zeros(length(PAT.Phase{kk}{jj}),nAmp);
        PAT.FlyMed.Mag          {kk,1}   	= zeros(length(PAT.Mag{kk}{jj}),nAmp);
      	PAT.FlyMed.Time        	{kk,1}   	= zeros(length(PAT.Time{kk}{jj}),nAmp);
        PAT.FlyMed.Pos          {kk,1}   	= zeros(length(PAT.Pos{kk}{jj}),nAmp);
        
        HEAD.FlyMed.Time        {kk,1}   	= zeros(length(HEAD.Time{kk}{jj}),nAmp);
        HEAD.FlyMed.Pos         {kk,1}   	= zeros(length(HEAD.Pos{kk}{jj}),nAmp);
        
     	HEAD.FlyMed.Freq        {kk,1}   	= zeros(length(HEAD.Freq{kk}{jj}),nAmp);
        HEAD.FlyMed.Phase       {kk,1}   	= zeros(length(HEAD.Phase{kk}{jj}),nAmp);
        HEAD.FlyMed.Mag         {kk,1}   	= zeros(length(HEAD.Mag{kk}{jj}),nAmp);
        
     	WING.FlyMed.Freq        {kk,1}   	= zeros(length(WING.Freq{kk}{jj}),nAmp);
        WING.FlyMed.Phase       {kk,1}   	= zeros(length(WING.Phase{kk}{jj}),nAmp);
        WING.FlyMed.Mag         {kk,1}   	= zeros(length(WING.Freq{kk}{jj}),nAmp);
        
        HEAD.FlyMed.Err.Freq    {kk,1}      = zeros(length(HEAD.GAIN{kk}{jj}),nAmp);
        HEAD.FlyMed.Err.Phase   {kk,1}      = zeros(length(HEAD.Phase{kk}{jj}),nAmp);
        HEAD.FlyMed.Err.Mag     {kk,1}      = zeros(length(HEAD.Mag{kk}{jj}),nAmp);
    end
end

for kk = 1:nFly
  	for jj = 1:nAmp
     	PAT.FlyMed.Freq         {kk}(:,jj) 	= median(PAT.Freq{kk}{jj},2);
        PAT.FlyMed.Phase        {kk}(:,jj) 	= median(PAT.Phase{kk}{jj},2);
        PAT.FlyMed.Mag          {kk}(:,jj) 	= median(PAT.Mag{kk}{jj},2);
        
     	PAT.FlyMed.Time         {kk}(:,jj) 	= median(PAT.Time{kk}{jj},2);
        PAT.FlyMed.Pos          {kk}(:,jj) 	= median(PAT.Pos{kk}{jj},2);
        
     	HEAD.FlyMed.Freq        {kk}(:,jj) 	= median(HEAD.Freq{kk}{jj},2);
        HEAD.FlyMed.Phase       {kk}(:,jj) 	= median(HEAD.Phase{kk}{jj},2);
        HEAD.FlyMed.Mag         {kk}(:,jj) 	= median(HEAD.Mag{kk}{jj},2);
        
    	WING.FlyMed.Freq        {kk}(:,jj) 	= median(WING.Freq{kk}{jj},2);
        WING.FlyMed.Phase       {kk}(:,jj)	= median(WING.Phase{kk}{jj},2);
        WING.FlyMed.Mag         {kk}(:,jj)	= median(WING.Mag{kk}{jj},2);
        
        HEAD.FlyMed.Err.Freq    {kk}(:,jj)  = median(HEAD.Err.Freq{kk}{jj},2);
        HEAD.FlyMed.Err.Phase   {kk}(:,jj)  = median(HEAD.Err.Phase{kk}{jj},2);
        HEAD.FlyMed.Err.Mag     {kk}(:,jj)  = median(HEAD.Err.Mag{kk}{jj},2);
        
        HEAD.FlyMed.Time        {kk}(:,jj) 	= median(HEAD.Time{kk}{jj},2);
        HEAD.FlyMed.Pos         {kk}(:,jj) 	= median(HEAD.Pos{kk}{jj},2);
    end
end

PAT.GrandMed.Freq  = median(cat(3,PAT.FlyMed.Freq{:}),3);
PAT.GrandMed.Phase = median(cat(3,PAT.FlyMed.Phase{:}),3);
PAT.GrandMed.Mag   = median(cat(3,PAT.FlyMed.Mag{:}),3);
PAT.GrandMed.Time  = median(cat(3,PAT.FlyMed.Time{:}),3);
PAT.GrandMed.Pos   = median(cat(3,PAT.FlyMed.Pos{:}),3);

HEAD.GrandMed.Freq  = median(cat(3,HEAD.FlyMed.Freq{:}),3);
HEAD.GrandMed.Phase = median(cat(3,HEAD.FlyMed.Phase{:}),3);
HEAD.GrandMed.Mag   = median(cat(3,HEAD.FlyMed.Mag{:}),3);

WING.GrandMed.Freq  = median(cat(3,WING.FlyMed.Freq{:}),3);
WING.GrandMed.Phase = median(cat(3,WING.FlyMed.Phase{:}),3);
WING.GrandMed.Mag   = median(cat(3,WING.FlyMed.Mag{:}),3);

HEAD.GrandMed.Err.Freq  = median(cat(3,HEAD.FlyMed.Err.Freq{:}),3);
HEAD.GrandMed.Err.Phase = median(cat(3,HEAD.FlyMed.Err.Phase{:}),3);
HEAD.GrandMed.Err.Mag   = median(cat(3,HEAD.FlyMed.Err.Mag{:}),3);

HEAD.GrandMed.Time  = median(cat(3,HEAD.FlyMed.Time{:}),3);
HEAD.GrandMed.Pos   = median(cat(3,HEAD.FlyMed.Pos{:}),3);
HEAD.GrandSTD.Pos   = std(cat(3,HEAD.FlyMed.Pos{:}),0,3);

%% Time Domain %%
%---------------------------------------------------------------------------------------------------------------------------------
clc
figure (76) ; clf
amp = 1:4;
pp = 1;
MR = cell(length(amp),1);
for kk = 1:nFly
    for jj = 1:length(amp)
        aa = uAmp(amp(jj));
        MD  = HEAD.GrandMed.Pos(:,amp(jj));
        SD  = HEAD.GrandSTD.Pos(:,amp(jj));
        DTA = HEAD.Pos{kk}{amp(jj)};

        for ii = 1:size(DTA,2)
            DT = DTA(:,ii);
            UP = MD + SD;
            LW = MD - SD;
            R = (DT<=UP) & (DT>=LW);
            pp = size(MR{jj},1) + 1;
            MR{jj}(pp,1) = mean(R);  % Value
            MR{jj}(pp,2) = kk; % Fly idx
            MR{jj}(pp,3) = aa; % Amplitude
         	MR{jj}(pp,4) = uFly(ii); % Trial (not actual)
            MR{jj}(pp,5) = uFly(kk); % Fly #
        end
        
        figure (76)
        subplot(length(amp),1,jj) ; hold on ; box on ; grid on ; title([num2str(uAmp(amp(jj))) '^{o}'])
%           	plot(HEAD.Time{kk}{amp(jj)},HEAD.Pos{kk}{amp(jj)},'Color',[0.5 0.5 0.5],'LineWidth',1)
            xlim([0 20])
            ylim([-20 20])
            xlabel('Time (s)')
            ylabel('deg')
    end 
end

for jj = 1:length(amp)
    figure (76)
    subplot(length(amp),1,jj)
    	plot(HEAD.GrandMed.Time(:,amp(jj)),PAT.GrandMed.Pos(:,amp(jj)),'b','LineWidth',3)
        g.patch = PatchSTD(HEAD.GrandMed.Pos(:,amp(jj)),HEAD.GrandSTD.Pos(:,amp(jj)),HEAD.GrandMed.Time(:,amp(jj)),...
        'k',[0.4 0.4 0.6]);
        plot(HEAD.GrandMed.Time(:,amp(jj)),HEAD.GrandMed.Pos(:,amp(jj)),'k','LineWidth',3)    
end

AMR = cell(length(MR),1);
for jj = 1:length(MR)
    D = MR{jj};
    [B,I]  = maxk(D(:,1),1);
    AMR{jj} = D(I,:);
    disp(AMR{jj})
    
    fi = AMR{jj}(:,2);
    aa = AMR{jj}(:,3);
    ai = amp(jj);
    ti = AMR{jj}(:,4);
    fn = AMR{jj}(:,5);
    
    for mm = 1:size(AMR{jj},1)
        figure (76)
        subplot(length(MR),1,jj)
            plot(HEAD.Time{fi(mm)}{ai}(:,ti(mm)),HEAD.Pos{fi(mm)}{ai}(:,ti(mm)),'LineWidth',2)
    end

end

%% Frequency Domain %%
%---------------------------------------------------------------------------------------------------------------------------------
% PAT FFT %
figure (46) ; clf
for kk = 1:nFly
    pp = 1;
    for jj = 1:nAmp
        figure (46)
            subplot(2,nAmp,pp) ; hold on
          	plot(PAT.FlyMed.Freq{kk}(:,jj),PAT.FlyMed.Mag{kk}(:,jj),'LineWidth',1)
            xlim([0.5 11.5])
%             ylim([0 1])
            box on ; grid on
            
            if kk==nFly
                plot(PAT.GrandMed.Freq(:,jj),PAT.GrandMed.Mag(:,jj),'k','LineWidth',3)
            end
            
            title(['Amp = ' num2str(uAmp(jj))])
            if pp==1 ; ylabel('PAT Mag') ; end

        figure (46)
            subplot(2,nAmp,pp+nAmp) ; hold on
          	plot(PAT.FlyMed.Freq{kk}(:,jj),PAT.FlyMed.Phase{kk}(:,jj),'LineWidth',1)
            xlim([0.5 11.5])
            ylim([-180 180])
            box on ; grid on
            
            if kk==nFly
                plot(PAT.GrandMed.Freq(:,jj),PAT.GrandMed.Phase(:,jj),'k','LineWidth',3)
            end
            
            xlabel('Frequency (Hz)')
         	if (pp+nAmp)==(1+nAmp) ; ylabel('PAT Phase (deg)') ; end

            pp = pp + 1;            
    end 
end

% HEAD FFT %
figure (47) ; clf
for kk = 1:nFly
    pp = 1;
    for jj = 1:nAmp
        figure (47)
            subplot(2,nAmp,pp) ; hold on
          	plot(HEAD.FlyMed.Freq{kk}(:,jj),HEAD.FlyMed.Mag{kk}(:,jj),'LineWidth',1)
            xlim([0.5 11.5])
%             ylim([0 1])
            box on ; grid on
            
            if kk==nFly
                plot(HEAD.GrandMed.Freq(:,jj),HEAD.GrandMed.Mag(:,jj),'k','LineWidth',3)
            end
            
            title(['Amp = ' num2str(uAmp(jj))])
            if pp==1 ; ylabel('Head Mag') ; end

        figure (47)
            subplot(2,nAmp,pp+nAmp) ; hold on
          	plot(HEAD.FlyMed.Freq{kk}(:,jj),HEAD.FlyMed.Phase{kk}(:,jj),'LineWidth',1)
            xlim([0.5 11.5])
            ylim([-180 180])
            box on ; grid on
            
            if kk==nFly
                plot(HEAD.GrandMed.Freq(:,jj),HEAD.GrandMed.Phase(:,jj),'k','LineWidth',3)
            end
            
            xlabel('Frequency (Hz)')
         	if (pp+nAmp)==(1+nAmp) ; ylabel('Head Phase (deg)') ; end

            pp = pp + 1;            
    end 
end

% WING (head-free) FFT %
figure (48) ; clf
for kk = 1:nFly
    pp = 1;
    for jj = 1:nAmp
        figure (48)
            subplot(2,nAmp,pp) ; hold on
          	plot(WING.FlyMed.Freq{kk}(:,jj),WING.FlyMed.Mag{kk}(:,jj),'LineWidth',1)
            xlim([0.5 11.5])
%             ylim([0 1])
            box on ; grid on
            
            if kk==nFly
                plot(WING.GrandMed.Freq(:,jj),WING.GrandMed.Mag(:,jj),'k','LineWidth',3)
            end
            
            title(['Amp = ' num2str(uAmp(jj))])
            if pp==1 ; ylabel('Wing Mag') ; end

        figure (48)
            subplot(2,nAmp,pp+nAmp) ; hold on
          	plot(WING.FlyMed.Freq{kk}(:,jj),WING.FlyMed.Phase{kk}(:,jj),'LineWidth',1)
            xlim([0.5 11.5])
            ylim([-180 180])
            box on ; grid on
            
            if kk==nFly
                plot(WING.GrandMed.Freq(:,jj),WING.GrandMed.Phase(:,jj),'k','LineWidth',3)
            end
            
            xlabel('Frequency (Hz)')
         	if (pp+nAmp)==(1+nAmp) ; ylabel('Wing Phase (deg)') ; end

            pp = pp + 1;            
    end 
end

% ERROR FFT %
figure (49) ; clf
for kk = 1:nFly
    pp = 1;
    for jj = 1:nAmp
        figure (49)
            subplot(2,nAmp,pp) ; hold on
          	plot(HEAD.FlyMed.Err.Freq{kk}(:,jj),HEAD.FlyMed.Err.Mag{kk}(:,jj),'LineWidth',1)
            xlim([0.5 11.5])
%             ylim([0 1])
            box on ; grid on
            
            if kk==nFly
                plot(HEAD.GrandMed.Err.Freq(:,jj),HEAD.GrandMed.Err.Mag(:,jj),'k','LineWidth',3)
            end
            
            title(['Amp = ' num2str(uAmp(jj))])
            if pp==1 ; ylabel(' Error Mag') ; end

        figure (49)
            subplot(2,nAmp,pp+nAmp) ; hold on
          	plot(HEAD.FlyMed.Err.Freq{kk}(:,jj),HEAD.FlyMed.Err.Phase{kk}(:,jj),'LineWidth',1)
            xlim([0.5 11.5])
            ylim([-180 180])
            box on ; grid on
            
            if kk==nFly
                plot(HEAD.GrandMed.Err.Freq(:,jj),HEAD.GrandMed.Err.Phase(:,jj),'k','LineWidth',3)
            end
            
            xlabel('Frequency (Hz)')
         	if (pp+nAmp)==(1+nAmp) ; ylabel('Error Phase (deg)') ; end

            pp = pp + 1;            
    end 
end

%% HEAD BODE %%
%---------------------------------------------------------------------------------------------------------------------------------
clc
if showplot.Bode
    figure (90) ; clf
end
for kk = 1:nFly
    for jj = 1:nAmp       
        HEAD.FlyMed.GAIN {kk,1}   = zeros(length(HEAD.GAIN{kk}{jj}),nAmp);
        HEAD.FlyMed.PHASE{kk,1}   = zeros(length(HEAD.PHASE{kk}{jj}),nAmp);
        HEAD.FlyMed.Freq {kk,1}   = zeros(length(HEAD.Freq{kk}{jj}),nAmp);
    end
end

for kk = 1:nFly
  	for jj = 1:nAmp       
        HEAD.FlyMed.GAIN{kk}(:,jj)  = median(HEAD.GAIN{kk}{jj},2);
        for ii = 1:size(HEAD.PHASE{kk}{jj},1)
            HEAD.FlyMed.PHASE{kk}(ii,jj) = circ_NewMean(HEAD.PHASE{kk}{jj}(ii,:));
        end
        HEAD.FlyMed.Freq{kk}(:,jj)  = median(HEAD.Freq{kk}{jj},2);
        
    	if showplot.Bode
        for ii = 1:size(HEAD.GAIN{kk}{jj},2)
            figure (90) % GAIN
                subIdx.gain = jj + 2*nAmp*(kk-1);
                subplot(2*nFly,nAmp,subIdx.gain) ; hold on
                plot(HEAD.Freq{kk}{jj}(:,ii),HEAD.GAIN{kk}{jj}(:,ii),'LineWidth',1)
                xlim([0.5 11.5])
                ylim([0 1])
                box on ; grid on
          	if subIdx.gain<=nAmp ; title(['Amp = ' num2str(uAmp(jj))]) ; end
            
            figure (90) % PHASE
                subIdx.phase = subIdx.gain + nAmp;
                subplot(2*nFly,nAmp,subIdx.phase) ; hold on
                plot(HEAD.Freq{kk}{jj}(:,ii),HEAD.PHASE{kk}{jj}(:,ii),'LineWidth',1)
                xlim([0.5 11.5])
                box on ; grid on
        end
        end
        
    	if showplot.Bode
        figure (90)
            subplot(2*nFly,nAmp,subIdx.gain) ; hold on
            plot(HEAD.FlyMed.Freq{kk}(:,jj),HEAD.FlyMed.GAIN{kk}(:,jj),'k','LineWidth',2)
        	if jj==1 ; ylabel(['Fly ' num2str(kk)]) ; end

        	subplot(2*nFly,nAmp,subIdx.phase) ; hold on
            plot(HEAD.FlyMed.Freq{kk}(:,jj),HEAD.FlyMed.PHASE{kk}(:,jj),'k','LineWidth',2)
        end
    end
end

HEAD.GrandMed.GAIN  = medfilt1(median(cat(3,HEAD.FlyMed.GAIN{:}),3),10);
HEAD.GrandMed.PHASE = medfilt1(rad2deg(median(cat(3,HEAD.FlyMed.PHASE{:}),3)),10);
HEAD.GrandMed.Freq  = median(cat(3,HEAD.FlyMed.Freq{:}),3);

HEAD.GrandSTD.GAIN  = std(cat(3,HEAD.FlyMed.GAIN{:}),0,3);
HEAD.GrandSTD.PHASE = rad2deg(std(cat(3,HEAD.FlyMed.PHASE{:}),0,3));

% BODE plot for all flies
figure (91) ; clf
for kk = 1:nFly
    pp = 1;
    for jj = 1:nAmp
        figure (91)
            subplot(2,nAmp,pp) ; hold on
          	plot(HEAD.FlyMed.Freq{kk}(:,jj),HEAD.FlyMed.GAIN{kk}(:,jj),'LineWidth',1)
            xlim([0.1 11.5])
            ylim([0 1])
            box on ; grid on
            
            if kk==nFly
                plot(HEAD.GrandMed.Freq(:,jj),HEAD.GrandMed.GAIN(:,jj),'k','LineWidth',3)
            end
            
            title([num2str(uAmp(jj)) ,'^{o}'])
            if pp==1 ; ylabel(' Head Gain') ; end

        figure (91)
            subplot(2,nAmp,pp+nAmp) ; hold on
          	plot(HEAD.FlyMed.Freq{kk}(:,jj),HEAD.FlyMed.PHASE{kk}(:,jj),'LineWidth',1)
            xlim([0.1 12])
%             ylim([-100 100])
            box on ; grid on
            
            if kk==nFly
                plot(HEAD.GrandMed.Freq(:,jj),HEAD.GrandMed.PHASE(:,jj),'k','LineWidth',3)
            end
            
            xlabel('Frequency (Hz)')
         	if (pp+nAmp)==(1+nAmp) ; ylabel('Head Phase (deg)') ; end

            pp = pp + 1;            
    end 
end

legList = cell(nFly,1);
for kk = 1:nFly
   legList{kk} = ['Fly ' num2str(kk)];  
end
legend(legList)

%% WING BODE minus ERROR %%
%---------------------------------------------------------------------------------------------------------------------------------
clc
for kk = 1:nFly
    for jj = 1:nAmp     	
        WING.FlyMed.Err.Gain {kk,1}  = zeros(length(WING.Err.GAIN{kk}{jj}),nAmp);
        WING.FlyMed.Err.Phase{kk,1}  = zeros(length(WING.Err.PHASE{kk}{jj}),nAmp);
%         WING.FlyMed.Err.Freq {kk,1}  = zeros(length(WING.Err.Freq{kk}{jj}),nAmp);
    end
end

for kk = 1:nFly
  	for jj = 1:nAmp        
        WING.FlyMed.Err.Gain{kk}(:,jj)  = median(WING.Err.GAIN{kk}{jj},2);
%         WING.FlyMed.Err.Phase{kk}(:,jj) = median(WING.Err.PHASE{kk}{jj},2);
        for ii = 1:size(WING.Err.PHASE{kk}{jj},1)
            WING.FlyMed.Err.Phase{kk}(ii,jj) = circ_NewMean(WING.Err.PHASE{kk}{jj}(ii,:));
        end
    end
end

WING.GrandMed.Err.Gain  = median(cat(3,WING.FlyMed.Err.Gain{:}),3);
WING.GrandMed.Err.Phase = medfilt1(rad2deg(median(cat(3,WING.FlyMed.Err.Phase{:}),3)),10);

WING.GrandSTD.Err.Gain  = std(cat(3,WING.FlyMed.Err.Gain{:}),0,3);
WING.GrandSTD.Err.Phase = rad2deg(std(cat(3,WING.FlyMed.Err.Phase{:}),0,3));

% BODE plot for all flies
figure (57) ; clf
for kk = 1:nFly
    pp = 1;
    for jj = 1:nAmp
        figure (57)
            subplot(2,nAmp,pp) ; hold on
          	plot(HEAD.FlyMed.Freq{kk}(:,jj),WING.FlyMed.Err.Gain{kk}(:,jj),'LineWidth',1)
            xlim([0.5 11.5])
            ylim([0 0.3])
            box on ; grid on
            
            if kk==nFly
                plot(HEAD.GrandMed.Freq(:,jj),WING.GrandMed.Err.Gain(:,jj),'k','LineWidth',3)
            end
            
            title(['Amp = ' num2str(uAmp(jj))])
            if pp==1 ; ylabel('Wing-Error Gain') ; end

        figure (57)
            subplot(2,nAmp,pp+nAmp) ; hold on
          	plot(HEAD.FlyMed.Freq{kk}(:,jj),WING.FlyMed.Err.Phase{kk}(:,jj),'LineWidth',1)
            xlim([0.5 11.5])
            ylim([-180 180])
            box on ; grid on
            
            if kk==nFly
                plot(HEAD.GrandMed.Freq(:,jj),WING.GrandMed.Err.Phase(:,jj),'k','LineWidth',3)
            end
            
            xlabel('Frequency (Hz)')
         	if (pp+nAmp)==(1+nAmp) ; ylabel('Wing-Error Phase (deg)') ; end

            pp = pp + 1;            
    end 
end

%% HEAD + WING %%
%---------------------------------------------------------------------------------------------------------------------------------
figure (400) ; clf ; 
set(gcf, 'Position',  [300, 200, 1200, 800])
set(gcf,'Color','w')
pp = 1;
amps = 2:4;
nPlot = length(amps);
for jj = amps
	subplot(2,nPlot,pp) ; hold on ; box on  ; title([num2str(uAmp(jj)) char(176)]) ; xlim([0.1 12]) ; grid on
    set(gca,'FontSize',28)
        yyaxis left ; set(gca,'YColor','b') ;if pp==1 ; ylabel(['Head Gain (' char(176) '/' char(176) ')']) ; end
            h.head.patch = PatchSTD(HEAD.GrandMed.GAIN(:,jj),HEAD.GrandSTD.GAIN(:,jj)./sqrt(nFly),HEAD.GrandMed.Freq(:,jj),...
                'b',[0.4 0.4 0.6]);
            h.head.med   = plot(HEAD.GrandMed.Freq(:,jj),HEAD.GrandMed.GAIN(:,jj),'-b','LineWidth',2);        
            ylim([0 1.0])
            if (pp==1)
                yticks([0 0.5 1])
            else
                yticklabels({})
            end
            xticklabels({})
            
        yyaxis right ;set(gca,'YColor','r') ; if pp==nPlot ; ylabel(['Wing Gain (V/' char(176) ')']) ; end
            h.wing.patch = PatchSTD(WING.GrandMed.Err.Gain(:,jj),WING.GrandSTD.Err.Gain(:,jj)./sqrt(nFly),HEAD.GrandMed.Freq(:,jj),...
                'r',[0.6 0.4 0.4]);
            h.wing.med = plot(HEAD.GrandMed.Freq(:,jj),WING.GrandMed.Err.Gain(:,jj),'-r','LineWidth',2);
            ylim([0 0.4])
            if (pp==3)
                yticks([0 0.2 0.4])
            else
                yticklabels({})
            end            
            if pp==nPlot ; legend([h.head.med,h.wing.med],'Head','Wings Head Free') ; end
    %-----------------------------------------------------------------------------------------------------------------------------
    subplot(2,nPlot,pp+nPlot) ; hold on ; box on ; grid on
    if (pp==2) ; xlabel('Frequency (Hz)') ; end
    set(gca,'FontSize',28)
        yyaxis left ; set(gca,'YColor','b') ; if pp==1 ; ylabel(['Head Phase (' char(176) ')']) ; end ; xlim([0.1 12])
            h.head.patch = PatchSTD(HEAD.GrandMed.PHASE(:,jj),HEAD.GrandSTD.PHASE(:,jj)./sqrt(nFly),HEAD.GrandMed.Freq(:,jj),...
                'b',[0.4 0.4 0.6]);
            h.head.med	 = plot(HEAD.GrandMed.Freq(:,jj),HEAD.GrandMed.PHASE(:,jj),'-b','LineWidth',2);        
            ylim([-90 90]) ; yticks([-90  0  90])
            if ~(pp==1) ; yticklabels({}) ;end
            
        yyaxis right ;set(gca,'YColor','r') ; if pp==nPlot ; ylabel(['Wing Phase (' char(176) ')']) ; end
            h.wing.patch = PatchSTD(WING.GrandMed.Err.Phase(:,jj),WING.GrandSTD.Err.Phase(:,jj)./sqrt(nFly),HEAD.GrandMed.Freq(:,jj),...
                'r',[0.6 0.4 0.4]);
            h.wing.med = plot(HEAD.GrandMed.Freq(:,jj),WING.GrandMed.Err.Phase(:,jj),'-r','LineWidth',2);
            ylim([-90 90]) ; yticks([-90  0  90])
            if ~(pp==3) ; yticklabels({}) ;end
            
    pp = pp + 1;
end

%% GAIN & PHASE %%
%---------------------------------------------------------------------------------------------------------------------------------
figure (401) ; clf ; set(gcf,'Color','w') ; set(gcf, 'Position',  [500, 100, 1100, 800])
figure (402) ; clf ; set(gcf,'Color','w') ; set(gcf, 'Position',  [400, 100, 1100, 800])
% set(gcf, 'Position',  [300, 200, 900, 450])
pp = 1;
amps = 2:4;
nPlot = length(amps);
for jj = amps
    figure (401)
	subplot(2,nPlot,pp) ; hold on ; box on  ; title([num2str(uAmp(jj)) char(176)]) ; xlim([0.1 12]) ; set(gca,'FontSize',28)
            if pp==1 ; ylabel(['Gain (' char(176) '/' char(176) ')']) ; end
            h.head.patch = PatchSTD(HEAD.GrandMed.GAIN(:,jj),HEAD.GrandSTD.GAIN(:,jj)./sqrt(nFly),HEAD.GrandMed.Freq(:,jj),...
                'b',[0.4 0.4 0.6]);
            h.head.med   = plot(HEAD.GrandMed.Freq(:,jj),HEAD.GrandMed.GAIN(:,jj),'-b','LineWidth',2);        
            ylim([0 1.0]) ; yticks([0 0.5 1])
            if ~(pp==1) ; yticklabels({}) ; end  
            xticklabels({})
            
	subplot(2,nPlot,pp+nPlot) ; hold on ; box on ; xlim([0.1 12]) ; set(gca,'FontSize',28)
          	if pp==1 ; ylabel(['Gain (V/' char(176) ')']) ; end       
            h.wing.patch = PatchSTD(WING.GrandMed.Err.Gain(:,jj),WING.GrandSTD.Err.Gain(:,jj)./sqrt(nFly),HEAD.GrandMed.Freq(:,jj),...
                'r',[0.6 0.4 0.4]);
            h.wing.med = plot(HEAD.GrandMed.Freq(:,jj),WING.GrandMed.Err.Gain(:,jj),'-r','LineWidth',2);
            ylim([0 0.3]) ; yticks([0 0.1 0.2 0.3])
            if ~(pp==1) ; yticklabels({}) ; end
            if  (pp==2) ; xlabel('Frequency (Hz)') ; end
            xticks([0.1 5 10])
            
            if pp==nPlot ; legend([h.head.med,h.wing.med],'Head','Wings Head Free') ; end
    %-----------------------------------------------------------------------------------------------------------------------------
    figure (402)
	subplot(2,nPlot,pp) ; hold on ; box on  ; title([num2str(uAmp(jj)) char(176)]) ; xlim([0.1 12]) ; set(gca,'FontSize',28)
            if pp==1 ; ylabel(['Phase (' char(176) ')']) ; end
            h.head.patch = PatchSTD(HEAD.GrandMed.PHASE(:,jj),HEAD.GrandSTD.PHASE(:,jj)./sqrt(nFly),HEAD.GrandMed.Freq(:,jj),...
                'b',[0.4 0.4 0.6]);
            h.head.med   = plot(HEAD.GrandMed.Freq(:,jj),HEAD.GrandMed.PHASE(:,jj),'-b','LineWidth',2);        
            ylim([-90 90]) ; yticks([-90 0 90])
            if ~(pp==1) ; yticklabels({}) ; end  
            xticklabels({})
            
	subplot(2,nPlot,pp+nPlot) ; hold on ; box on  ; xlim([0.1 12]) ; set(gca,'FontSize',28)
          	if pp==1 ; ylabel(['Phase (' char(176) ')']) ; end       
            h.wing.patch = PatchSTD(WING.GrandMed.Err.Phase(:,jj),WING.GrandSTD.Err.Phase(:,jj)./sqrt(nFly),HEAD.GrandMed.Freq(:,jj),...
                'r',[0.6 0.4 0.4]);
            h.wing.med = plot(HEAD.GrandMed.Freq(:,jj),WING.GrandMed.Err.Phase(:,jj),'-r','LineWidth',2);
            ylim([-90 90]) ; yticks([-90 0 90])
            if ~(pp==1) ; yticklabels({}) ; end
            if  (pp==2) ; xlabel('Frequency (Hz)') ; end
            xticks([0.1 5 10])
            
        	if pp==nPlot ; legend([h.head.med,h.wing.med],'Head','Wings Head Free') ; end 
    pp = pp + 1;
end

%% LOG PLOT GAIN & PHASE %%
%---------------------------------------------------------------------------------------------------------------------------------
figure (403) ; clf ; set(gcf,'Color','w') ; set(gcf, 'Position',  [500, 100, 1100, 800])
figure (404) ; clf ; set(gcf,'Color','w') ; set(gcf, 'Position',  [300, 100, 1100, 800])
% set(gcf, 'Position',  [300, 200, 900, 450])
clc
pp = 1;
amps = 2:4;
nPlot = length(amps);
for jj = amps
    figure (403)
	subplot(2,nPlot,pp) ; semilogx(0) ;  hold on ; box on  ; title([num2str(uAmp(jj)) char(176)]) ; set(gca,'FontSize',28)
            if pp==1 ; ylabel(['Head Gain (' char(176) '/' char(176) ')']) ; end
            
            DD = HEAD.GrandMed.GAIN(:,jj);
            SD = HEAD.GrandSTD.GAIN(:,jj);
            FD = HEAD.GrandMed.Freq(:,jj);
            
            h.head.patch = PatchSTD(DD(2:end),SD(2:end)./sqrt(nFly),FD(2:end),...
                'b',[0.4 0.4 0.6]);
            h.head.med   = semilogx(HEAD.GrandMed.Freq(:,jj),HEAD.GrandMed.GAIN(:,jj),'-b','LineWidth',2);        
            ylim([0 1.0]) ; yticks([0 0.5 1])
            
            if ~(pp==1) ; yticklabels({}) ; end  
            xticklabels({})
            
	subplot(2,nPlot,pp+nPlot) ; semilogx(0); hold on ; box on ; set(gca,'FontSize',28)
          	if pp==1 ; ylabel(['Wing Gain (V/' char(176) ')']) ; end       
            
          	DD = WING.GrandMed.Err.Gain(:,jj);
            SD = WING.GrandSTD.Err.Gain(:,jj);
            FD = HEAD.GrandMed.Freq(:,jj);
            
            h.wing.patch = PatchSTD(DD(2:end),SD(2:end)./sqrt(nFly),FD(2:end),...
                'r',[0.6 0.4 0.4]);
            h.wing.med = semilogx(HEAD.GrandMed.Freq(:,jj),WING.GrandMed.Err.Gain(:,jj),'-r','LineWidth',2);
            
            ylim([0 0.3]) ; yticks([0 0.1 0.2 0.3])
            if ~(pp==1) ; yticklabels({}) ; end
            if  (pp==2) ; xlabel('Frequency (Hz)') ; end
            
            if pp==nPlot ; legend([h.head.med,h.wing.med],'Head','Wings Head Free') ; end
    %-----------------------------------------------------------------------------------------------------------------------------
    figure (404)
	subplot(2,nPlot,pp) ; semilogx(0) ; hold on ; box on  ; title([num2str(uAmp(jj)) char(176)]) ; set(gca,'FontSize',28)
            if pp==1 ; ylabel(['Head Phase (' char(176) ')']) ; end
            
          	DD = HEAD.GrandMed.PHASE(:,jj);
            SD = HEAD.GrandSTD.PHASE(:,jj);
            FD = HEAD.GrandMed.Freq(:,jj);
            
            h.head.patch = PatchSTD(DD(2:end),SD(2:end)./sqrt(nFly),FD(2:end),...
                'b',[0.4 0.4 0.6]);
            h.head.med   = semilogx(HEAD.GrandMed.Freq(:,jj),HEAD.GrandMed.PHASE(:,jj),'-b','LineWidth',2);   
                   
            ylim([-90 90]) ; yticks([-90 0 90])
            if ~(pp==1) ; yticklabels({}) ; end  
            xticklabels({})
            
	subplot(2,nPlot,pp+nPlot) ; semilogx(0) ; hold on ; box on  ; set(gca,'FontSize',28)
          	if pp==1 ; ylabel(['Wing Phase (' char(176) ')']) ; end
          	
            DD = WING.GrandMed.Err.Phase(:,jj);
            SD = WING.GrandSTD.Err.Phase(:,jj);
            FD = HEAD.GrandMed.Freq(:,jj);
            
            h.wing.patch = PatchSTD(DD(2:end),SD(2:end)./sqrt(nFly),FD(2:end),...
                'r',[0.6 0.4 0.4]);
            h.wing.med = semilogx(HEAD.GrandMed.Freq(:,jj), WING.GrandMed.Err.Phase(:,jj),'-r','LineWidth',2);
            
            ylim([-90 90]) ; yticks([-90 0 90])
            if ~(pp==1) ; yticklabels({}) ; end
            if  (pp==2) ; xlabel('Frequency (Hz)') ; end
            
    pp = pp + 1;
end


%% SYSTEM ID %%
%---------------------------------------------------------------------------------------------------------------------------------
clc
% figure (28) ; clf
figure (27) ; clf
SI = 10;
EI = 211;
for jj = 1:1
    FreqData    = HEAD.GrandMed.Freq (SI:EI,jj);
    GainData    = HEAD.GrandMed.GAIN (SI:EI,jj);
    PhaseData   = HEAD.GrandMed.PHASE(SI:EI,jj);
    
    zfr = GainData.*exp(1i*PhaseData*pi/180);
    Ts = 1/200;
    gfr = idfrd(zfr,FreqData,Ts);
    
    figure (27)
    bode(gfr)
    
    m1(jj) = oe(gfr,[2 2 1]) % Discrete-time Output error (transfer function) model
    ms(jj) = ssest(gfr) % Continuous-time state-space model with default choice of order
    mproc(jj) = procest(gfr,'P2UDZ') % 2nd-order, continuous-time model with underdamped poles
    compare(gfr,m1(jj),ms(jj),mproc(jj))
    
%     figure (28)
%     subplot(2,1,1) ; hold on ; title('Gain') ; box on
%         plot(FreqData,GainData,'LineWidth',2)
%     
% 	subplot(2,1,2) ; hold on ; title('Phase') ; box on
%         plot(FreqData,PhaseData,'LineWidth',2)
end




%% HEAD COHERENCE %%
%---------------------------------------------------------------------------------------------------------------------------------
clc
if showplot.Coher
    figure (70) ; clf
end
HEAD.FlyMean.COHR.Mag = [];
HEAD.FlyMean.COHR.Freq = [];
for kk = 1:nFly
    for jj = 1:nAmp
        HEAD.FlyMean.COHR.Mag {kk,1}  = zeros(length(HEAD.COHR.Mag{kk}{jj}),nAmp);
      	HEAD.FlyMean.COHR.Freq {kk,1}  = zeros(length(HEAD.COHR.Freq{kk}{jj}),nAmp);
    end
end

for kk = 1:nFly
  	for jj = 1:nAmp
        HEAD.FlyMean.COHR.Mag  {kk}(:,jj) = mean(HEAD.COHR.Mag{kk}{jj},2);
        HEAD.FlyMean.COHR.Freq {kk}(:,jj) = mean(HEAD.COHR.Freq{kk}{jj},2);
        if showplot.Coher
        for ii = 1:size(HEAD.GAIN{kk}{jj},2)
            figure (70) % GAIN
                subIdx.gain = jj + nAmp*(kk-1);
                subplot(nFly,nAmp,subIdx.gain) ; hold on
                plot(HEAD.COHR.Freq{kk}{jj}(:,ii),HEAD.COHR.Mag{kk}{jj}(:,ii),'LineWidth',1)
                xlim([0.5 11.5])
                box on ; grid on
          	if subIdx.gain<=nAmp ; title(['Amp = ' num2str(uAmp(jj))]) ; end
        end
        end
        if showplot.Coher
        figure (70)
            subplot(nFly,nAmp,subIdx.gain) ; hold on
            plot(HEAD.FlyMean.COHR.Freq{kk}(:,jj),HEAD.FlyMean.COHR.Mag{kk}(:,jj),'k','LineWidth',2)
        	if jj==1 ; ylabel(['Fly ' num2str(kk)]) ; end
        end
    end
end

HEAD.GrandMean.COHR.Mag  = mean(cat(3,HEAD.FlyMean.COHR.Mag{:}),3);
HEAD.GrandMean.COHR.Freq = mean(cat(3,HEAD.FlyMean.COHR.Freq{:}),3);

% COHERENCE for all flies
figure (71) ; clf
for kk = 1:nFly
    pp = 1;
    for jj = 1:nAmp
        figure (71)
            subplot(1,nAmp,pp) ; hold on
          	plot(HEAD.FlyMean.COHR.Freq{kk}(:,jj),HEAD.FlyMean.COHR.Mag{kk}(:,jj),'LineWidth',1)
            xlim([0.5 11.5])
            ylim([0 1])
            box on ; grid on
            
            if kk==nFly
                plot(HEAD.GrandMean.COHR.Freq(:,jj),HEAD.GrandMean.COHR.Mag(:,jj),'k','LineWidth',3)
            end
            
            title(['Amp = ' num2str(uAmp(jj))])
            if pp==1 ; ylabel(' Head Coherence') ; end

            pp = pp + 1;            
    end 
end

%% WING COHERENCE %%
%---------------------------------------------------------------------------------------------------------------------------------
clc
if showplot.Coher
    figure (60) ; clf
end
WING.FlyMean.COHR.Mag = [];
WING.FlyMean.COHR.Freq = [];
for kk = 1:nFly
    for jj = 1:nAmp
        WING.FlyMean.COHR.Mag {kk,1}  = zeros(length(WING.COHR.Mag{kk}{jj}),nAmp);
      	WING.FlyMean.COHR.Freq {kk,1}  = zeros(length(WING.COHR.Freq{kk}{jj}),nAmp);
    end
end

for kk = 1:nFly
  	for jj = 1:nAmp
        WING.FlyMean.COHR.Mag  {kk}(:,jj) = mean(WING.COHR.Mag{kk}{jj},2);
        WING.FlyMean.COHR.Freq {kk}(:,jj) = mean(WING.COHR.Freq{kk}{jj},2);
        if showplot.Coher
        for ii = 1:size(WING.GAIN{kk}{jj},2)
            figure (60) % GAIN
                subIdx.gain = jj + nAmp*(kk-1);
                subplot(nFly,nAmp,subIdx.gain) ; hold on
                plot(WING.COHR.Freq{kk}{jj}(:,ii),WING.COHR.Mag{kk}{jj}(:,ii),'LineWidth',1)
                xlim([0.5 11.5])
                box on ; grid on
          	if subIdx.gain<=nAmp ; title(['Amp = ' num2str(uAmp(jj))]) ; end
        end
        end
        if showplot.Coher
        figure (60)
            subplot(nFly,nAmp,subIdx.gain) ; hold on
            plot(WING.FlyMean.COHR.Freq{kk}(:,jj),WING.FlyMean.COHR.Mag{kk}(:,jj),'k','LineWidth',2)
        	if jj==1 ; ylabel(['Fly ' num2str(kk)]) ; end
        end
    end
end

WING.GrandMean.COHR.Mag  = mean(cat(3,WING.FlyMean.COHR.Mag{:}),3);
WING.GrandMean.COHR.Freq = mean(cat(3,WING.FlyMean.COHR.Freq{:}),3);

% COHERENCE for all flies
figure (61) ; clf
for kk = 1:nFly
    pp = 1;
    for jj = 1:nAmp
        figure (61)
            subplot(1,nAmp,pp) ; hold on
          	plot(WING.FlyMean.COHR.Freq{kk}(:,jj),WING.FlyMean.COHR.Mag{kk}(:,jj),'LineWidth',1)
            xlim([0.5 11.5])
            ylim([0 1])
            box on ; grid on
            
            if kk==nFly
                plot(WING.GrandMean.COHR.Freq(:,jj),WING.GrandMean.COHR.Mag(:,jj),'k','LineWidth',3)
            end
            
            title(['Amp = ' num2str(uAmp(jj))])
            if pp==1 ; ylabel(' Wing Coherence') ; end

            pp = pp + 1;            
    end 
end
end
%% FUNCTION:	fitPanel
function [xi,ti] = fitPanel(xPos,time,time_new)
%---------------------------------------------------------------------------------------------------------------------------------
% fitPanel: fits curve to panel panel pattern position
    % INPUTS:
        % xPos: pattern position in deg
        % t_p: time vector
 	% OUTPUTS
        % ti: new time vector
        % xi: fitted data
        % timeMean: interrpolated time vector
        % xPosMean: interrpolated data points
%---------------------------------------------------------------------------------------------------------------------------------      
    % Find when panel position changes (transition)
    pp = 1;
    for jj = 2:length(xPos)
        dP = xPos(jj) - xPos(jj-1);
        if abs(dP) > 1
            markerIndex(pp) = jj-1;
            markerPos(pp) = xPos(jj-1);
            markerTime(pp) = time(jj-1);
            pp = pp + 1;
        end
    end
    
    % Find the mean values between transisiton points
    indexMean = nan(1,length(markerIndex));
    timeMean = nan(1,length(markerIndex));
    xPosMean = nan(1,length(markerIndex));
    for jj = 2:length(markerIndex)
        indexMean(jj-1) = round(mean(markerIndex(jj-1):1:markerIndex(jj)));
        timeMean(jj-1)  = time   (indexMean(jj-1));
        xPosMean(jj-1)  = xPos  (indexMean(jj-1));
    end

    timeMean = [time(1) timeMean(1:end-1) time_new(end)]; % time when data is between two transisiton points
    xPosMean = [xPos(1) xPosMean(1:end-1) xPos(end)]; % pattern position when data is between two transisiton points
    
    ti = timeMean;
    xi = xPosMean;
    
    xi = spline(timeMean',xPosMean',time_new); % spline interrpolated data
    ti = time_new;

% 	figure (1) ; clf ; hold on
%         plot(markerTime,markerPos,'g*')
%         plot(time,xPos,'b')
%         plot(timeMean,xPosMean,'r*-');
%         plot(time_new,xi,'-k')
end
%---------------------------------------------------------------------------------------------------------------------------------
%% FUNCTION:    FFT
function [Fv, Mag , Phase] = FFT(t,x)
%---------------------------------------------------------------------------------------------------------------------------------
% FFT: Transforms time domian data to frequency domain
    % INPUTS:
        % t: time data
        % x: pos/vel data
	% OUTPUTS
        % Fv: frequency vector
        % Mag: magniude of fft
        % Phase: phase of fft
%---------------------------------------------------------------------------------------------------------------------------------      
    Fs = 1/(mean(diff(t)));                 % Sampling frequency [Hz]
    L = length(t);                          % Length of signal
    Fn = Fs/2;                           	% Nyquist Frequency
    fts = fft(x)/L;                        	% Normalised Fourier Transform
    Fv = (linspace(0, 1, fix(L/2)+1)*Fn)';  % Frequency Vector
    Iv = 1:length(Fv);                  	% Index Vector
    
    Mag = abs(fts(Iv))*2;                   % Magnitude
    Phase = (angle(fts(Iv)));      	% Phase
%---------------------------------------------------------------------------------------------------------------------------------
end
%---------------------------------------------------------------------------------------------------------------------------------
%% FUNCTION:    PatchSTD
function [h] = PatchSTD(dataMean,dataSTD,dataX,edgeColor,faceColor)
%---------------------------------------------------------------------------------------------------------------------------------
% PatchSTD: Plots patch for data mean
    % INPUTS:
        % dataMean: mean of data
        % dataSTD: std of data
      	% dataX: x-axis points
        % edgeColor: patch edge color
        % faceColor: pathc area color
	% OUTPUTS
        % h.patch: patch handle
%---------------------------------------------------------------------------------------------------------------------------------      
    uE = dataMean + 1*dataSTD;
    lE = dataMean - 1*dataSTD;
    yP = [lE;flipud(uE)];
    xP = [dataX;flipud(dataX)];
    h = patch(xP,yP,1,'facecolor',faceColor,'edgecolor',edgeColor);
    alpha(h,0.5)
%---------------------------------------------------------------------------------------------------------------------------------
end