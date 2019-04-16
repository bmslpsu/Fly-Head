showplot.Time  = 0;
showplot.Freq  = 0;

savename = 'Chirp_Walking_Data';

rootdir = 'W:\Research\Walking Chirp mat';
[FILES, PATH] = uigetfile({'*.mat', 'DAQ-files'}, 'Select files', rootdir, 'MultiSelect','on');
FILES = cellstr(FILES)'; % if only one character array >> store in cell

[Original,Index,Number,Unique] = GetFileData(FILES,'Fly','Trial','Amp');

%% Get Data %%
%---------------------------------------------------------------------------------------------------------------------------------
TIME                = cell(Number{1,1},1);

HEAD.Pos            = cell(Number{1,1},1);
HEAD.Freq           = cell(Number{1,1},1);
HEAD.Mag            = cell(Number{1,1},1);
HEAD.Phase          = cell(Number{1,1},1);

PAT.Pos             = cell(Number{1,1},1);
PAT.Freq            = cell(Number{1,1},1);
PAT.Mag             = cell(Number{1,1},1);
PAT.Phase           = cell(Number{1,1},1);

BODE.Mag.HeadPat    = cell(Number{1,1},1);
BODE.Freq           = cell(Number{1,1},1);
BODE.Phase.HeadPat  = cell(Number{1,1},1);

for kk = 1:Number{1,1} % fly
    for jj = 1:Number{1,3} % Amplitude
            FLY.Time{kk,1}{jj,1}            = [];
            
            HEAD.Pos{kk,1}{jj,1}            = [];
            HEAD.Freq{kk,1}{jj,1}           = [];
            HEAD.Mag{kk,1}{jj,1}            = [];
            HEAD.Phase{kk,1}{jj,1}          = [];
            HEAD.COHR.Freq{kk,1}{jj,1}      = [];
            HEAD.COHR.Mag{kk,1}{jj,1}       = [];
            
            PAT.Pos{kk,1}{jj,1}             = [];
            PAT.Freq{kk,1}{jj,1}            = [];
            PAT.Mag{kk,1}{jj,1}             = [];
            PAT.Phase{kk,1}{jj,1}           = [];
            
            BODE.Mag.HeadPat{kk,1}{jj,1}    = [];
            BODE.Freq{kk,1}{jj,1}           = [];
            BODE.Phase.HeadPat{kk,1}{jj,1}  = [];
            
            
    end
end
fly = [];
FlyState = [];
AI = [];
fly.Fc = 20;
span = 20:1:2100;
for kk = 1:Number{1,end} % all trials
    filename = fullfile(PATH,FILES{kk}); % full file name
    load(filename,'FlyState','AI','VidTime') % load in fly kinematics & arena voltages
    
    vid.time        = VidTime;
    fly.Time        = FlyState{:,1};
    fly.Ts          = mean(diff(fly.Time));
    fly.Fs          = 1/fly.Ts; 
	[b,a]           = butter(2,fly.Fc/(fly.Fs/2),'low'); % 2nd-order low-pass butterworth filter
    head.Pos        = filtfilt(b,a,FlyState{:,2});
	head.Vel        = [diff(head.Pos)./fly.Ts ; 0];
    pat.Time        = AI{:,1};
    pat.Ts          = mean(diff(pat.Time));
    pat.Fs          = 1/pat.Ts;
 	pat.Pos         = AI{:,2};

 	head.Pos        = interp1(fly.Time, head.Pos , vid.time, 'nearest'); % interpolate pattern y-pos to match fly
 	head.Vel        = interp1(fly.Time, head.Vel , vid.time, 'nearest'); % interpolate pattern y-pos to match fly
 	pat.Pos         = interp1(pat.Time, pat.Pos  , vid.time, 'nearest'); % interpolate pattern x-pos to match fly
 	
    fly.Time        = fly.Time  (span);
    head.Pos        = head.Pos  (span);
	head.Vel        = head.Vel  (span);
 	pat.Pos         = pat.Pos   (span);
    
    pat.Pos         = 3.75.*round(pat.Pos.*96/5);       % from voltage to degrees
    

    [head.Freq , head.Mag, head.Phase]   = FFT(fly.Time,head.Pos);
    [pat.Freq  , pat.Mag  , pat.Phase]   = FFT(fly.Time,pat.Pos);
    SI = 1; EI = 238;
    head.Freq   = head.Freq(SI:EI);
    head.Mag    = head.Mag(SI:EI);
    head.Phase  = head.Phase(SI:EI);
    pat.Freq    = pat.Freq(SI:EI);
    pat.Mag     = pat.Mag(SI:EI);
    pat.Phase   = pat.Phase(SI:EI);
    
    [head.cohr.mag, head.cohr.f] = mscohere(pat.Pos ,head.Pos,[],[],head.Freq,fly.Fs);  % Calculate coherence
	
    
    head.Pos    = rad2deg(head.Pos);        % head angle to degrees
    
    FLY.Time{Index{kk,1}}{Index{kk,3},1}(:,end+1)               = fly.Time;
    
    HEAD.Pos{Index{kk,1}}{Index{kk,3},1}(:,end+1)               = head.Pos;
    HEAD.Freq{Index{kk,1}}{Index{kk,3},1}(:,end+1)              = head.Freq;
    HEAD.Mag{Index{kk,1}}{Index{kk,3},1}(:,end+1)               = head.Mag;
    HEAD.Phase{Index{kk,1}}{Index{kk,3},1}(:,end+1)             = head.Phase;
    
    HEAD.COHR.Freq  {Index{kk,1}}{Index{kk,3},1}(:,end+1)       = head.cohr.f;
    HEAD.COHR.Mag   {Index{kk,1}}{Index{kk,3},1}(:,end+1)       = head.cohr.mag;
  	
    PAT.Pos{Index{kk,1}}{Index{kk,3},1}(:,end+1)                = pat.Pos ;
    PAT.Freq{Index{kk,1}}{Index{kk,3},1}(:,end+1)               = pat.Freq ;
    PAT.Mag{Index{kk,1}}{Index{kk,3},1}(:,end+1)                = pat.Mag ;
    PAT.Phase{Index{kk,1}}{Index{kk,3},1}(:,end+1)              = pat.Phase ;
    
    BODE.Mag.HeadPat{Index{kk,1}}{Index{kk,3},1}(:,end+1)       = medfilt1(head.Mag./pat.Mag,1);
    BODE.Freq{Index{kk,1}}{Index{kk,3},1}(:,end+1)              = pat.Freq;
    BODE.Phase.HeadPat{Index{kk,1}}{Index{kk,3},1}(:,end+1)     = -1*medfilt1(pat.Phase - head.Phase,1);
   
  %--------------------------------------------------------------------------------------------------------------------------------- 
% Plot previews

colmn = 4; 
    if showplot.Time
        figure (100)
            subplot(ceil(Number{1,4}/colmn),colmn,kk) ; hold on
                title(['Fly ' num2str(Original{kk,1})])
                plot(fly.Time,pat.Pos,'k')
                plot(fly.Time,head.Pos,'b')
                box on
    end      
    if showplot.Freq
        figure(104)
            subplot(ceil(Number{1,4}/colmn),colmn,kk) ; hold on
                title(['Fly ' num2str(Original{kk,1})])
                plot(pat.Freq,  pat.Mag ,'k')
                plot(head.Freq, head.Mag,'b')
                xlim([0.5 11.5])
                box on
                hold off
        figure (105)
            subplot(ceil(Number{1,4}/colmn),colmn,kk) ; hold on
                title(['Fly ' num2str(Original{kk,1})])
                plot(pat.Freq,  pat.Phase ,'k')
                plot(head.Freq, head.Phase ,'b')
                xlim([0.5 11.5])
                box on    
    end
end



%% Averages

clc
% for kk = 1:Number{1,1}
%     for jj = 1:Number{1,3}
%         
%         TIME.FlyMed         	   	= zeros(length(TIME{kk}{jj}(:,1)),Number{1,3});
%         
%      	PAT.FlyMed.Freq         {kk,1}   	= zeros(length(PAT.Freq{kk}{jj}(:,1)),Number{1,3});
%         PAT.FlyMed.Phase        {kk,1}   	= zeros(length(PAT.Phase{kk}{jj(:,1)}),Number{1,3});
%         PAT.FlyMed.Mag          {kk,1}   	= zeros(length(PAT.Mag{kk}{jj}(:,1)),Number{1,3});
%         PAT.FlyMed.Pos          {kk,1}   	= zeros(length(PAT.Pos{kk}{jj}(:,1)),Number{1,3});
%         
%      	HEAD.FlyMed.Freq        {kk,1}   	= zeros(length(HEAD.Freq{kk}{jj}(:,1)),Number{1,3});
%         HEAD.FlyMed.Phase       {kk,1}   	= zeros(length(HEAD.Phase{kk}{jj}(:,1)),Number{1,3});
%         HEAD.FlyMed.Mag         {kk,1}   	= zeros(length(HEAD.Mag{kk}{jj}(:,1)),Number{1,3});
%         HEAD.FlyMed.Pos         {kk,1}   	= zeros(length(HEAD.Pos{kk}{jj}(:,1)),Number{1,3});
%         
%         BODE.FlyMed.Freq        {kk,1}      = zeros(length(HEAD.Freq{kk}{jj}(:,1)),Number{1,3});
%         BODE.FlyMed.Phase       {kk,1}      = zeros(length(HEAD.Phase{kk}{jj}(:,1)),Number{1,3});
%         BODE.FlyMed.Mag         {kk,1}      = zeros(length(HEAD.Mag{kk}{jj}(:,1)),Number{1,3});
%     end
% end


for kk = 1:Number{1,1}
  	for jj = 1:Number{1,3}
        
        FLY.Time.FlyMed         {kk}{:,jj}   	  = median(FLY.Time{kk}{jj},2);
                
     	PAT.FlyMed.Freq         {kk}(:,jj)        = median(PAT.Freq{kk}{jj},2);
        PAT.FlyMed.Phase        {kk}(:,jj)        = median(PAT.Phase{kk}{jj},2);
        PAT.FlyMed.Mag          {kk}(:,jj)        = median(PAT.Mag{kk}{jj},2);
        PAT.FlyMed.Pos          {kk}(:,jj)        = median(PAT.Pos{kk}{jj},2);
        
     	HEAD.FlyMed.Freq        {kk}(:,jj)        = median(HEAD.Freq{kk}{jj},2);
        HEAD.FlyMed.Phase       {kk}(:,jj)        = median(HEAD.Phase{kk}{jj},2);
        HEAD.FlyMed.Mag         {kk}(:,jj)        = median(HEAD.Mag{kk}{jj},2);
        HEAD.FlyMed.Pos         {kk}(:,jj)        = median(HEAD.Pos{kk}{jj},2);
        
        BODE.FlyMed.Freq        {kk,1}(:,jj)      = median(BODE.Freq{kk}{jj},2);
        BODE.FlyMed.Phase       {kk,1}(:,jj)      = median(BODE.Phase{kk}{jj},2);
        BODE.FlyMed.Mag         {kk,1}(:,jj)      = median(BODE.Mag{kk}{jj}, 2);
    end
end

FLY.Time.GrandMed      = median(cat(3,FLY.Time.FlyMed{1}), 3);
PAT.GrandMed.Freq  = median(cat(3,PAT.FlyMed.Freq{:}),3);
PAT.GrandMed.Phase = median(cat(3,PAT.FlyMed.Phase{:}),3);
PAT.GrandMed.Mag   = median(cat(3,PAT.FlyMed.Mag{:}),3);
PAT.GrandMed.Pos   = median(cat(3,PAT.FlyMed.Pos{:}),3);

HEAD.GrandMed.Freq  = median(cat(3,HEAD.FlyMed.Freq{:}),3);
HEAD.GrandMed.Phase = median(cat(3,HEAD.FlyMed.Phase{:}),3);
HEAD.GrandMed.Mag   = median(cat(3,HEAD.FlyMed.Mag{:}),3);
HEAD.GrandMed.Pos   = median(cat(3,HEAD.FlyMed.Pos{:}),3);
HEAD.GrandSTD.Pos   = std(cat(3,HEAD.FlyMed.Pos{:}),0,3);

BODE.GrandMed.Freq  = median(cat(3,BODE.FlyMed.Freq{:}),3);
BODE.GrandMed.Phase = median(cat(3,BODE.FlyMed.Phase{:}),3);
BODE.GrandMed.Mag   = median(cat(3,BODE.FlyMed.Mag{:}),3);


