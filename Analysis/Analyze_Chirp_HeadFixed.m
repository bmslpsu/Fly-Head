function [] = Analyze_Chirp_HeadFixed()
% Analyze_Chirp_HeadFixed: Summary of this function goes here
%   Detailed explanation goes here
%   INPUTS:
%       root: root directory
%   OUTPUTS:
%
%---------------------------------------------------------------------------------------------------------------------------------
showplot.Time = 0;
showplot.Freq = 1;
showplot.Bode = 0;
showplot.Coher = 0;
%% Setup Directories %%
%---------------------------------------------------------------------------------------------------------------------------------
root.pat = 'E:\Experiment_HeadExcitation\Chirp\HeadFixed\';

% Select files
[FILES, PATH.pat] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root.pat, 'MultiSelect','on');
FILES = FILES';

%% Process File Data %%
%---------------------------------------------------------------------------------------------------------------------------------
clearvars -except FILES PATH root showplot
% close all
clc
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
clear PAT WING head pat wings Amplitude
% close all
SI = 1;
EI = 241;
Amplitude = zeros(nTrial,1);
% Preallocate data cells %
for kk = 1:nFly % # of flys
    for jj = 1:nAmp % # of amplitudes      
        PAT.Time{kk,1}{jj,1}        = [];
        PAT.Pos{kk,1}{jj,1}         = [];
        PAT.Freq{kk,1}{jj,1}        = [];
        PAT.Mag{kk,1}{jj,1}          = [];
        PAT.Phase{kk,1}{jj,1}       = [];
        PAT.Vel{kk,1}{jj,1}         = [];
        
        WING.Time{kk,1}{jj,1}       = [];
        WING.Pos{kk,1}{jj,1}        = [];
        WING.Freq{kk,1}{jj,1}       = [];
        WING.Mag{kk,1}{jj,1}         = [];
        WING.Phase{kk,1}{jj,1}      = [];
        WING.GAIN{kk,1}{jj,1}       = [];
        WING.PHASE{kk,1}{jj,1}      = [];        
        WING.COHR.Freq{kk,1}{jj,1}  = [];
        WING.COHR.Mag{kk,1}{jj,1}   = [];
    end
end
% Save time & frequency domain data in cells %
for kk = 1:nTrial
    clear pat data t_v t_p  % clear temp variables
    %-----------------------------------------------------------------------------------------------------------------------------
    % Load DAQ data %
	load([PATH.pat   FILES{kk}]); % load pattern x-position
	%-----------------------------------------------------------------------------------------------------------------------------
    % Get wing data from DAQ %
  	K = 1; % WBA factor
    wing.TimeFull  = t_p; % wing time
    
  	% Setup filter for wings %
    wing.Fs = 1/mean(diff(wing.TimeFull)); % sampling frequency [Hz]
    wing.Fc = 20; % cutoff frequency [Hz]
    [b,a] = butter(2,wing.Fc/(wing.Fs/2)); % butterworth filter
    
	wing.Left       = filtfilt(b,a,(data(:,4))); % left wing
    wing.Right      = filtfilt(b,a,(data(:,5))); % right wing
    wing.PosFull    = (K*(wing.Left - wing.Right)); % WBA (L-R)
    wing.f          = 100*(data(:,6)); % wing beat frequency
    
 	if min(wing.f)<180
        fprintf('Low WBF: Fly %i Trial %i \n',Fly(kk),Trial(kk))
        continue
    end
    
    wing.PosFull = wing.PosFull - mean(wing.PosFull); % subtract mean
    wing.Pos  = wing.PosFull;
    wing.Time = t_p;
    wing.Vel = [diff(wing.Pos)./diff(wing.Time) ; 0];
	%-----------------------------------------------------------------------------------------------------------------------------
    % Get pattern data from DAQ %
    pat.TimeFull = t_p;
    pat.PosFull  = 3.75*round((96/10)*(data(:,2)' - mean(data(:,2))))'; % get pattern x-pos: subtract mean and convert to deg
    pat.Pos  = pat.PosFull;
	pat.Time = pat.TimeFull;
	pat.Vel = [diff(pat.Pos)./diff(pat.Time) ; 0];
    Amplitude(kk,1) = max(pat.Pos);
    %-----------------------------------------------------------------------------------------------------------------------------
    % Convert head, wings, & pattern data to frequency domain %
    [wing.Freq , wing.Mag , wing.Phase]  = FFT(wing.Time,wing.Pos);
    [pat.Freq  , pat.Mag  , pat.Phase ]  = FFT(pat.Time,pat.Pos);
    
 	wing.Freq   = wing.Freq(SI:EI);
    wing.Mag     = wing.Mag(SI:EI);
    wing.Phase  = wing.Phase(SI:EI);
	
    pat.Freq   = pat.Freq(SI:EI);
    pat.Mag     = pat.Mag(SI:EI);
    pat.Phase  = pat.Phase(SI:EI);   
    %-----------------------------------------------------------------------------------------------------------------------------
    % Calculate BODE gain % phase %       
	wing.bode.mag = wing.Mag./pat.Mag; % wing gain
    wing.bode.phase = medfilt1(-(pat.Phase - wing.Phase),20); % wing phase
    %-----------------------------------------------------------------------------------------------------------------------------
    % Calculate coherence %
	[wing.cohr.mag,wing.cohr.f] = mscohere(pat.Pos ,wing.Pos ,[],[],wing.Freq ,wing.Fs);
    %-----------------------------------------------------------------------------------------------------------------------------
    % Store data in cells %  
    % Pattern
	PAT.Time  	{idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = pat.Time;
	PAT.Pos    	{idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = pat.Pos;
	PAT.Vel    	{idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = pat.Vel;

    PAT.Freq  	{idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = pat.Freq;
	PAT.Mag      {idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = pat.Mag;
	PAT.Phase	{idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = pat.Phase;
    
    % Wings
	WING.Pos	{idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = wing.Pos;
	WING.Time	{idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = wing.Time;

	WING.Freq  {idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = wing.Freq;
	WING.Mag    {idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = wing.Mag;
	WING.Phase	{idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = wing.Phase;
    
	WING.GAIN   {idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = wing.bode.mag;
	WING.PHASE  {idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = wing.bode.phase;
    
	WING.COHR.Freq  {idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = wing.cohr.f;
    WING.COHR.Mag   {idxFly(kk),1}{idxAmp(kk),1}(:,end+1) = wing.cohr.mag;    
    %-----------------------------------------------------------------------------------------------------------------------------
    colmn = 4;
    if showplot.Time
        figure (100)
            subplot(ceil(nTrial/colmn),colmn,kk) ; hold on
                title(['Fly ' num2str(Fly(kk)) ' Amp ' num2str(Amp(kk))])
                plot(pat.Time,pat.Pos,'k')
                plot(wing.Time,wing.Pos,'b')
                box on
    end
    if showplot.Freq
        figure (104)
            subplot(ceil(nTrial/colmn),colmn,kk) ; hold on
                title(['Fly ' num2str(Fly(kk)) ' Amp ' num2str(Amp(kk))])
                plot(pat.Freq,  pat.Mag ,'g')
                plot(wing.Freq, wing.Mag ,'b')
                xlim([0.5 11.5])
                box on
        figure (105)
            subplot(ceil(nTrial/colmn),colmn,kk) ; hold on
                title(['Fly ' num2str(Fly(kk)) ' Amp ' num2str(Amp(kk))])
                plot(pat.Freq,  pat.Phase ,'g')
                plot(wing.Freq, wing.Phase ,'b')
                xlim([0.5 11.5])
                box on    
    end
    
    if showplot.Bode
        figure (109)
            subplot(ceil(nTrial/colmn),colmn,kk) ; hold on
                title(['Fly ' num2str(Fly(kk)) ' Amp ' num2str(Amp(kk))])
                plot(wing.Freq, wing.bode.mag ,'k')
                xlim([0.5 11.5])
                box on
        figure (110)
            subplot(ceil(nTrial/colmn),colmn,kk) ; hold on
                title(['Fly ' num2str(Fly(kk)) ' Amp ' num2str(Amp(kk))])
                plot(wing.Freq, wing.bode.phase ,'k')
                xlim([0.5 11.5])
                box on
    end
    %-----------------------------------------------------------------------------------------------------------------------------
end
disp('DONE')
%% Averages %%
%---------------------------------------------------------------------------------------------------------------------------------
for kk = 1:nFly
    for jj = 1:nAmp
        PAT.FlyMed.Time  	{kk,1} 	= zeros(length(PAT.Time{kk}{jj}),nAmp);
        PAT.FlyMed.Pos    	{kk,1} 	= zeros(length(PAT.Pos{kk}{jj}),nAmp);
     	PAT.FlyMed.Freq    	{kk,1}	= zeros(length(PAT.Freq{kk}{jj}),nAmp);
        PAT.FlyMed.Phase  	{kk,1} 	= zeros(length(PAT.Phase{kk}{jj}),nAmp);
        PAT.FlyMed.Mag    	{kk,1} 	= zeros(length(PAT.Mag{kk}{jj}),nAmp);
        
        WING.FlyMed.Time 	{kk,1} 	= zeros(length(WING.Time{kk}{jj}),nAmp);
        WING.FlyMed.Pos  	{kk,1}	= zeros(length(WING.Pos{kk}{jj}),nAmp);
     	WING.FlyMed.Freq 	{kk,1} 	= zeros(length(WING.Freq{kk}{jj}),nAmp);
        WING.FlyMed.Phase	{kk,1} 	= zeros(length(WING.Phase{kk}{jj}),nAmp);
        WING.FlyMed.Mag   	{kk,1}	= zeros(length(WING.Freq{kk}{jj}),nAmp);
        
        WING.FlyMed.GAIN   {kk,1}  = zeros(length(WING.GAIN{kk}{jj}),nAmp);
        WING.FlyMed.PHASE  {kk,1}  = zeros(length(WING.PHASE{kk}{jj}),nAmp);
    end
end

for kk = 1:nFly
  	for jj = 1:nAmp
        PAT.FlyMed.Time  	{kk}(:,jj) 	= median(PAT.Time{kk}{jj},2);
        PAT.FlyMed.Pos    	{kk}(:,jj) 	= median(PAT.Pos{kk}{jj},2);
     	PAT.FlyMed.Freq  	{kk}(:,jj) 	= median(PAT.Freq{kk}{jj},2);
        PAT.FlyMed.Phase  	{kk}(:,jj) 	= median(PAT.Phase{kk}{jj},2);
        PAT.FlyMed.Mag   	{kk}(:,jj) 	= median(PAT.Mag{kk}{jj},2);

        WING.FlyMed.Time   	{kk}(:,jj) 	= median(WING.Time{kk}{jj},2);
        WING.FlyMed.Pos  	{kk}(:,jj) 	= median(WING.Pos{kk}{jj},2);
    	WING.FlyMed.Freq    {kk}(:,jj) 	= median(WING.Freq{kk}{jj},2);
        WING.FlyMed.Phase	{kk}(:,jj)	= median(WING.Phase{kk}{jj},2);
        WING.FlyMed.Mag    	{kk}(:,jj)	= median(WING.Mag{kk}{jj},2);
    end
end

PAT.GrandMed.Time  = median(cat(3,PAT.FlyMed.Time{:}),3);
PAT.GrandMed.Pos   = median(cat(3,PAT.FlyMed.Pos{:}),3);
PAT.GrandMed.Freq  = median(cat(3,PAT.FlyMed.Freq{:}),3);
PAT.GrandMed.Phase = median(cat(3,PAT.FlyMed.Phase{:}),3);
PAT.GrandMed.Mag   = median(cat(3,PAT.FlyMed.Mag{:}),3);

WING.GrandMed.Time  = median(cat(3,WING.FlyMed.Time{:}),3);
WING.GrandMed.Pos   = median(cat(3,WING.FlyMed.Pos{:}),3);
WING.GrandMed.Freq  = median(cat(3,WING.FlyMed.Freq{:}),3);
WING.GrandMed.Phase = median(cat(3,WING.FlyMed.Phase{:}),3);
WING.GrandMed.Mag   = median(cat(3,WING.FlyMed.Mag{:}),3);
WING.GrandSTD.Pos   = std(cat(3,WING.FlyMed.Pos{:}),0,3);

%% WING TIME DOMAIN
%---------------------------------------------------------------------------------------------------------------------------------
figure (101) ; clf
amp = 1:4;
set(gcf,'Position',[100 100 1100 800/4])
% suptitle('Wings: Head Fixed')
for kk = 1:nFly
  	for jj = 1:nAmp
    	subplot(nAmp,1,jj) ; hold on ; title(['Amp = ' num2str(uAmp(jj))])
        h = plot(WING.Time{kk}{jj},WING.Pos{kk}{jj},'Color',[0.5 0.5 0.5]);
        for ii = 1:length(h)
            h(ii).Color(4) = 0.1;
        end
        h = plot(WING.FlyMed.Time{kk}(:,jj),WING.FlyMed.Pos{kk}(:,jj),'Color',[0.5 0.5 0.5]);
        h.Color(4) = 0.5;
    end
end
clear h
for jj = 1:nAmp
	subplot(nAmp,1,jj) ; hold on ; title([num2str(uAmp(amp(jj))) '$^{\circ}$'],'Interpreter','latex','FontSize',15)
%       	plot(PAT.GrandMed.Time(:,jj), PAT.GrandMed.Pos(:,jj),'b','LineWidth',1);
        h.med = plot(WING.GrandMed.Time(:,jj),WING.GrandMed.Pos(:,jj),'k','LineWidth',3);
        h.patch = PatchSTD(WING.GrandMed.Pos(:,amp(jj)),WING.GrandSTD.Pos(:,amp(jj)),WING.GrandMed.Time(:,amp(jj)),...
        'k',[0.4 0.4 0.6]);
    
     	xlim([0 20])
        ylim([-4 4])
%         yticks([-20:10:20])
        if jj==nAmp
            xlabel('Time (s)','Interpreter','latex','FontSize',15)
        end
        if jj~=nAmp
            xticks([0])
            xticklabels([''])
        end
        ylabel('$\Delta$ WBA$(V)$','Interpreter','latex','FontSize',15)
end

%% WING BODE %%
%---------------------------------------------------------------------------------------------------------------------------------
for kk = 1:nFly
  	for jj = 1:nAmp
        WING.FlyMean.Gain{kk}(:,jj)  = nanmean(WING.GAIN{kk}{jj},2);
        WING.FlyMean.Phase{kk}(:,jj) = nanmean(WING.PHASE{kk}{jj},2);
        for ii = 1:size(WING.PHASE{kk}{jj},1)
            WING.FlyMed.Phase{kk}(ii,jj) = circ_NewMean(WING.PHASE{kk}{jj}(ii,:));
        end
        WING.FlyMean.Freq{kk}(:,jj)  = nanmean(WING.Freq{kk}{jj},2);
        
        WING.FlyMed.Gain{kk}(:,jj)  = nanmedian(WING.GAIN{kk}{jj},2);
        WING.FlyMed.Phase{kk}(:,jj) = nanmedian(WING.PHASE{kk}{jj},2);
        WING.FlyMed.Freq{kk}(:,jj)  = nanmedian(WING.Freq{kk}{jj},2);
        
        if showplot.Bode
        for ii = 1:size(WING.GAIN{kk}{jj},2)
            figure (80) % GAIN
                subIdx.gain = jj + 2*nAmp*(kk-1);
                subplot(2*nFly,nAmp,subIdx.gain) ; hold on
                plot(WING.Freq{kk}{jj}(:,ii),WING.GAIN{kk}{jj}(:,ii),'LineWidth',1)
                xlim([0.5 11.5])
                box on ; grid on
          	if subIdx.gain<=nAmp ; title(['Amp = ' num2str(uAmp(jj))]) ; end
            
            figure (80) % PHASE
                subIdx.phase = subIdx.gain + nAmp;
                subplot(2*nFly,nAmp,subIdx.phase) ; hold on
                plot(WING.Freq{kk}{jj}(:,ii),WING.PHASE{kk}{jj}(:,ii),'LineWidth',1)
                xlim([0.5 11.5])
                box on ; grid on
        end
        end
        if showplot.Bode
        figure (80)
            subplot(2*nFly,nAmp,subIdx.gain) ; hold on
            plot(WING.FlyMean.Freq{kk}(:,jj),WING.FlyMean.Gain{kk}(:,jj),'k','LineWidth',2)
        	if jj==1 ; ylabel(['Fly ' num2str(kk)]) ; end

        	subplot(2*nFly,nAmp,subIdx.phase) ; hold on
            plot(WING.FlyMean.Freq{kk}(:,jj),WING.FlyMean.Phase{kk}(:,jj),'k','LineWidth',2)
        end
    end
end

WING.GrandMean.Gain  = nanmean(cat(3,WING.FlyMean.Gain{:}),3);
WING.GrandMean.Phase = nanmean(cat(3,WING.FlyMean.Phase{:}),3);
WING.GrandMean.Freq  = nanmean(cat(3,WING.FlyMean.Freq{:}),3);

WING.GrandMed.Gain  = medfilt1(median(cat(3,WING.FlyMed.Gain{:}),3),10);
WING.GrandMed.Phase = rad2deg(median(cat(3,WING.FlyMed.Phase{:}),3));
WING.GrandMed.Freq  = median(cat(3,WING.FlyMed.Freq{:}),3);

WING.GrandSTD.Gain  = std(cat(3,WING.FlyMean.Gain{:}),0,3);
WING.GrandSTD.Phase = rad2deg(std(cat(3,WING.FlyMean.Phase{:}),0,3));

% BODE plot for all flies
figure (81) ; clf
for kk = 1:nFly
    pp = 1;
    for jj = 1:nAmp
        figure (81)
            subplot(2,nAmp,pp) ; hold on
          	plot(WING.FlyMean.Freq{kk}(:,jj),WING.FlyMean.Gain{kk}(:,jj),'LineWidth',1)
            xlim([0.5 11.5])
            box on ; grid on
            
            if kk==nFly
                plot(WING.GrandMean.Freq(:,jj),WING.GrandMean.Gain(:,jj),'k','LineWidth',3)
            end
            
            title(['Amp = ' num2str(uAmp(jj))])
            if pp==1 ; ylabel('Wing Gain') ; end

        figure (81)
            subplot(2,nAmp,pp+nAmp) ; hold on
          	plot(WING.FlyMean.Freq{kk}(:,jj),WING.FlyMean.Phase{kk}(:,jj),'LineWidth',1)
            xlim([0.5 11.5])
            box on ; grid on
            
            if kk==nFly
                plot(WING.GrandMean.Freq(:,jj),WING.GrandMean.Phase(:,jj),'k','LineWidth',3)
            end
            
            xlabel('Frequency (Hz)')
         	if (pp+nAmp)==(1+nAmp) ; ylabel('Wing Phase (deg)') ; end

            pp = pp + 1;            
    end 
end

%% GAIN & PHASE %%
%---------------------------------------------------------------------------------------------------------------------------------
pp = 1;
amps = 2:4;
nPlot = length(amps);
for jj = amps
    figure (401)
	subplot(2,nPlot,pp+nPlot) ; hold on 
            h.wing.patch = PatchSTD(WING.GrandMed.Gain(:,jj),WING.GrandSTD.Gain(:,jj)./sqrt(nFly),WING.GrandMed.Freq(:,jj),...
                'k',[0.8 0.8 0.8]);
            h.wing.med = plot(WING.GrandMed.Freq(:,jj),WING.GrandMed.Gain(:,jj),'-k','LineWidth',2);
            
    %-----------------------------------------------------------------------------------------------------------------------------
    figure (402)
	subplot(2,nPlot,pp+nPlot) ; hold on
            h.wing.patch = PatchSTD(WING.GrandMed.Phase(:,jj),WING.GrandSTD.Phase(:,jj)./sqrt(nFly),WING.GrandMed.Freq(:,jj),...
                'k',[0.8 0.8 0.8]);
            h.wing.med = plot(WING.GrandMed.Freq(:,jj),WING.GrandMed.Phase(:,jj),'-k','LineWidth',2);

    pp = pp + 1;
end

%% WING COHERENCE %%
%---------------------------------------------------------------------------------------------------------------------------------
clc
figure (60) ; clf
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
        WING.FlyMean.COHR.Mag  {kk}(:,jj) = nanmean(WING.COHR.Mag{kk}{jj},2);
        WING.FlyMean.COHR.Freq {kk}(:,jj) = nanmean(WING.COHR.Freq{kk}{jj},2);

        for ii = 1:size(WING.GAIN{kk}{jj},2)
            figure (60) % GAIN
                subIdx.gain = jj + nAmp*(kk-1);
                subplot(nFly,nAmp,subIdx.gain) ; hold on
                plot(WING.COHR.Freq{kk}{jj}(:,ii),WING.COHR.Mag{kk}{jj}(:,ii),'LineWidth',1)
                xlim([0.5 11.5])
                box on ; grid on
          	if subIdx.gain<=nAmp ; title(['Amp = ' num2str(uAmp(jj))]) ; end

        end
               
        figure (60)
            subplot(nFly,nAmp,subIdx.gain) ; hold on
            plot(WING.FlyMean.COHR.Freq{kk}(:,jj),WING.FlyMean.COHR.Mag{kk}(:,jj),'k','LineWidth',2)
        	if jj==1 ; ylabel(['Fly ' num2str(kk)]) ; end
    end
end

WING.GrandMean.COHR.Mag  = nanmean(cat(3,WING.FlyMean.COHR.Mag{:}),3);
WING.GrandMean.COHR.Freq = nanmean(cat(3,WING.FlyMean.COHR.Freq{:}),3);

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
    Phase = (angle(fts(Iv)));               % Phase
%---------------------------------------------------------------------------------------------------------------------------------
end
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