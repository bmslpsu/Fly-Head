showplot.Time  = 1;
showplot.Freq  = 1;
npeaks = 6 ;             % Number of sine waves in pattern
%% Setup Directories %%
%---------------------------------------------------------------------------------------------------------------------------------
root.pat = '/Users/akb427/Documents/MATLAB/Research/Sum-of-Sine';
root.head = [root.pat '/Vid/Angles/'];

% Select files
[FILES, PATH.head] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root.head, 'MultiSelect','on');
FILES = FILES';

PATH.pat = uigetdir(root.pat);
PATH.pat = [PATH.pat '/'];

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

fprintf('Total Flies: %i \n',nFly) ; fprintf('Total Trials: %i \n',nTrial)
T = cell2table(trialFly,'VariableNames',{'Fly','Trials'});
disp(T)

%% Get Data %%
%---------------------------------------------------------------------------------------------------------------------------------
disp('Loading Data...')
clear HEAD PAT WING head pat wings 
SI = 1;
EI = 241;
% Preallocate data cells %
for kk = 1:nFly % # of flys
        HEAD.Time{kk,1}       = [];
        HEAD.Pos{kk,1}        = [];
        HEAD.Freq{kk,1}       = [];
        HEAD.Mag{kk,1}     	  = [];
        HEAD.Phase{kk,1}      = [];
        HEAD.Peak.Mag{kk,1}   = [];
        HEAD.Peak.Freq{kk,1}  = [];
        HEAD.Peak.Phase{kk,1} = [];
        
        ERROR.Time{kk,1}        = [];
        ERROR.Pos{kk,1}         = [];
        ERROR.Freq{kk,1}        = [];
        ERROR.Mag{kk,1}         = [];
        ERROR.Phase{kk,1}       = [];
        ERROR.Peak.Mag{kk,1}    = [];
        ERROR.Peak.Freq{kk,1}   = [];
        ERROR.Peak.Phase{kk,1}  = [];
        
        PAT.Time{kk,1}       = [];
        PAT.Pos{kk,1}        = [];
        PAT.Freq{kk,1}       = [];
        PAT.Mag{kk,1}        = [];
        PAT.Phase{kk,1}      = [];
        PAT.Peak.Mag{kk,1}   = [];
        PAT.Peak.Freq{kk,1}  = [];
        PAT.Peak.Phase{kk,1} = [];
        
        WING.Time{kk,1}      = [];
        WING.Pos{kk,1}       = [];
        WING.Freq{kk,1}      = [];
        WING.Mag{kk,1}    	 = [];
        WING.Phase{kk,1}     = [];
        WING.Peak.Mag{kk,1}         =[];
        WING.Peak.Freq{kk,1}        =[];
        WING.Peak.Phase{kk,1}       =[];
        
        BODE.Mag.HeadPat{kk,1}          =[];
        BODE.Mag.WingError{kk,1}        =[];
        BODE.Freq{kk,1}                 =[];
        BODE.Phase.HeadPat{kk,1}        =[];
        BODE.Phase.WingError{kk,1}      =[];
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
    
    interval = floor(length(wing.Full.Time)/length(head.Time));
    wing.Pos  = wing.Full.Pos(1:interval:end-1);
    wing.Time = wing.Full.Time(1:interval:end-1);
           
 	wing.f = 100*(data(:,6)); % wing beat frequency
    if min(wing.f)<180
        fprintf('Low WBF: Fly %i Trial %i \n',Fly(kk),Trial(kk))
    end
    
%-----------------------------------------------------------------------------------------------------------------------------
% Get pattern data from DAQ %
    pat.Full.Time = t_p;
    pat.Full.Pos  = 3.75*round((96/10)*(data(:,2)' - mean(data(:,2))))'; % get pattern x-pos: subtract mean and convert to deg

    interval = floor(length(pat.Full.Time)/length(head.Time));
    pat.Pos  = pat.Full.Pos(1:interval:end-1);
	pat.Time = pat.Full.Time(1:interval:end-1);
    
%-----------------------------------------------------------------------------------------------------------------------------
      % Compute Error in Time Domain %
      error.Pos = pat.Pos-head.Pos;
      error.Time = pat.Time;
%-----------------------------------------------------------------------------------------------------------------------------
     % Convert head, wings, & pattern data to frequency domain %
    [head.Freq , head.Mag, head.Phase]   = FFT(head.Time,head.Pos);
    [wing.Freq , wing.Mag , wing.Phase]  = FFT(wing.Time,wing.Pos);
    [pat.Freq  , pat.Mag  , pat.Phase ]  = FFT(pat.Time,pat.Pos);
    [error.Freq, error.Mag, error.Phase] = FFT(error.Time,error.Pos);
    
    head.Freq   = head.Freq(SI:EI);
    head.Mag    = head.Mag(SI:EI);
    head.Phase  = head.Phase(SI:EI);
    
    wing.Freq   = wing.Freq(SI:EI);
    wing.Mag    = wing.Mag(SI:EI);
    wing.Phase  = wing.Phase(SI:EI);
	
    pat.Freq   = pat.Freq(SI:EI);
    pat.Mag    = pat.Mag(SI:EI);
    pat.Phase  = pat.Phase(SI:EI);
    
    error.Freq   = error.Freq(SI:EI);
    error.Mag    = error.Mag(SI:EI);
    error.Phase  = error.Phase(SI:EI);
%-----------------------------------------------------------------------------------------------------------------------------
    % Find Peak Frequencies & Phases %

    [pat.Peak.Mag,pat.Peak.Loc] = maxk(pat.Mag,npeaks);
    [pat.Peak.Loc,ind1] = sort(pat.Peak.Loc);
    pat.Peak.Mag = pat.Peak.Mag(ind1);
    pat.Peak.Freq = pat.Freq(pat.Peak.Loc);
    pat.Peak.Phase = pat.Phase(pat.Peak.Loc);
    
    [head.Peak.Freq, head.Peak.Mag, head.Peak.Phase] = SOSPeakFinder(npeaks, pat.Peak.Freq, head.Freq, head.Mag, head.Phase );
    [error.Peak.Freq, error.Peak.Mag, error.Peak.Phase] = SOSPeakFinder(npeaks, pat.Peak.Freq, error.Freq, error.Mag, error.Phase);
    [wing.Peak.Freq, wing.Peak.Mag, wing.Peak.Phase] = SOSPeakFinder(npeaks, pat.Peak.Freq, wing.Freq, wing.Mag, wing.Phase);
    
%-----------------------------------------------------------------------------------------------------------------------------
    % Store data in cells %
    % Head
	HEAD.Time     {idxFly(kk),1}(:,end+1)       = head.Time;
	HEAD.Pos  	  {idxFly(kk),1}(:,end+1)       = head.Pos;
	HEAD.Freq   {idxFly(kk),1}(:,end+1)         = head.Freq;
	HEAD.Mag 	{idxFly(kk),1}(:,end+1)         = head.Mag;
	HEAD.Phase	{idxFly(kk),1}(:,end+1)         = head.Phase;
    HEAD.Peak.Freq {idxFly(kk),1}(:,end+1)      = head.Peak.Freq;
    HEAD.Peak.Mag {idxFly(kk),1}(:,end+1)       = head.Peak.Mag;
    HEAD.Peak.Phase {idxFly(kk),1}(:,end+1)     = head.Peak.Phase;
    
    % Pattern
	PAT.Time  	{idxFly(kk),1}(:,end+1)     = pat.Time;
	PAT.Pos    	{idxFly(kk),1}(:,end+1)     = pat.Pos;
    PAT.Freq  	{idxFly(kk),1}(:,end+1)     = pat.Freq;
	PAT.Mag    	{idxFly(kk),1}(:,end+1)     = pat.Mag;
	PAT.Phase	{idxFly(kk),1}(:,end+1)     = pat.Phase;
    PAT.Peak.Freq {idxFly(kk),1}(:,end+1)   = pat.Peak.Freq;
    PAT.Peak.Mag {idxFly(kk),1}(:,end+1)    = pat.Peak.Mag;
    PAT.Peak.Phase {idxFly(kk),1}(:,end+1)  = pat.Peak.Phase;
    
    % Wings
	WING.Pos {idxFly(kk),1}(:,end+1)        = wing.Pos;
	WING.Time {idxFly(kk),1}(:,end+1)       = wing.Time;
	WING.Freq {idxFly(kk),1}(:,end+1)       = wing.Freq;
	WING.Mag {idxFly(kk),1}(:,end+1)        = wing.Mag;
	WING.Phase {idxFly(kk),1}(:,end+1)      = wing.Phase;
    HEAD.Peak.Freq {idxFly(kk),1}(:,end+1)  = wing.Peak.Freq;
    HEAD.Peak.Mag {idxFly(kk),1}(:,end+1)   = wing.Peak.Mag;
    HEAD.Peak.Phase {idxFly(kk),1}(:,end+1) = wing.Peak.Phase;
    
    % Error
    ERROR.Time{idxFly(kk),1}(:,end+1)       = error.Time;
    ERROR.Pos{idxFly(kk),1}(:,end+1)        = error.Pos;
    ERROR.Freq{idxFly(kk),1}(:,end+1)       = error.Freq;
    ERROR.Mag{idxFly(kk),1}(:,end+1)        = error.Mag;
    ERROR.Phase{idxFly(kk),1}(:,end+1)      = error.Phase;
    ERROR.Peak.Mag{idxFly(kk),1}(:,end+1)   = error.Peak.Mag;
    ERROR.Peak.Freq{idxFly(kk),1}(:,end+1)  = error.Peak.Freq;
    ERROR.Peak.Phase{idxFly(kk),1}(:,end+1) = error.Peak.Phase;
     
%-----------------------------------------------------------------------------------------------------------------------------
    colmn = 4; 
    if showplot.Time
        figure (100)
            subplot(ceil(nTrial/colmn),colmn,kk) ; hold on
                title(['Fly ' num2str(Fly(kk))])
                plot(pat.Time,pat.Pos,'k')
                plot(head.Time,head.Pos,'b')
                plot(wing.Time, wing.Pos, 'r')
                legend('pat','head','wing')
                box on
    end
               
    if showplot.Freq
        figure(104)
            subplot(ceil(nTrial/colmn),colmn,kk) ; hold on
                title(['Fly ' num2str(Fly(kk))])
                plot(pat.Freq,  pat.Mag ,'k')
                plot(pat.Peak.Freq, pat.Peak.Mag, 'ko')
                plot(head.Freq, head.Mag,'b')
                plot(head.Peak.Freq, head.Peak.Mag, 'bo')
             	plot(wing.Freq, wing.Mag,'r')
                plot(wing.Peak.Freq, wing.Peak.Mag, 'ro')
                xlim([0.5 11.5])
                box on
                hold off
        figure (105)
            subplot(ceil(nTrial/colmn),colmn,kk) ; hold on
                title(['Fly ' num2str(Fly(kk))])
                plot(pat.Freq,  pat.Phase ,'k')
                plot(pat.Peak.Freq, pat.Peak.Phase, 'ko')
                plot(head.Freq, head.Phase ,'b')
                plot(head.Peak.Freq, head.Peak.Mag, 'bo')
                plot(wing.Freq, wing.Phase,'r')
                plot(wing.Peak.Freq, wing.Peak.Mag, 'ro')
                xlim([0.5 11.5])
                box on    
    end
%-----------------------------------------------------------------------------------------------------------------------------
    % Calculate Bode Data %
    BODE.Mag.HeadPat{idxFly(kk),1}(:,end+1)             = head.Peak.Mag./pat.Peak.Mag;
    BODE.Mag.WingError{idxFly(kk),1}(:,end+1)           = wing.Peak.Mag./error.Peak.Mag;
    BODE.Freq{idxFly(kk),1}(:,end+1)                    = pat.Peak.Freq;
    BODE.Phase.HeadPat{idxFly(kk),1}(:,end+1)           = head.Peak.Phase- pat.Peak.Phase;
    BODE.Phase.WingError{idxFly(kk),1}(:,end+1)         = wing.Peak.Phase - error.Peak.Phase;
%-----------------------------------------------------------------------------------------------------------------------------
end


disp('DONE')

    