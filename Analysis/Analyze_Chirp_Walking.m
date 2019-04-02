showplot.Time  = 1;
showplot.Freq  = 1;

rootdir = 'C:\Users\amb7595\Documents\Experimental Data\Walking Chirp';
[FILES, PATH] = uigetfile({'*.mat', 'DAQ-files'}, 'Select files', rootdir, 'MultiSelect','on');
FILES = cellstr(FILES)'; % if only one character array >> store in cell

[~,Index,Number,Unique] = GetFileData(FILES,'Fly','Trial','Amp');

%% Get Data %%
%---------------------------------------------------------------------------------------------------------------------------------
TIME        = cell(Number{1,1},1);

HEAD.Pos    = cell(Number{1,1},1);
HEAD.Freq   = cell(Number{1,1},1);
HEAD.Mag    = cell(Number{1,1},1);
HEAD.Phase  = cell(Number{1,1},1);

PAT.Pos  	= cell(Number{1,1},1);
PAT.Freq    = cell(Number{1,1},1);
PAT.Mag     = cell(Number{1,1},1);
PAT.Phase   = cell(Number{1,1},1);

BODE.Mag.HeadPat   = cell(Number{1,1},1);
BODE.Freq          = cell(Number{1,1},1);
BODE.Phase.HeadPat = cell(Number{1,1},1);

for kk = 1:Number{1,1} % fly
    for jj = 1:Number{1,3} % Amplitude
            TIME{kk,1}{jj,1}       = [];
            HEAD.Pos{kk,1}{jj,1}   = [];
            HEAD.Freq{kk,1}{jj,1}  = [];
            HEAD.Mag{kk,1}{jj,1}   = [];
            HEAD.Phase{kk,1}{jj,1} = [];
            
            PAT.Pos{kk,1}{jj,1}   = [];
            PAT.Freq{kk,1}{jj,1}  = [];
            PAT.Mag{kk,1}{jj,1}   = [];
            PAT.Phase{kk,1}{jj,1} = [];
            
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
    
    vid.time            = VidTime;
    fly.time            = FlyState{:,1};
    fly.Ts              = mean(diff(fly.time));
    fly.Fs              = 1/fly.Ts; 
	[b,a]               = butter(2,fly.Fc/(fly.Fs/2),'low'); % 2nd-order low-pass butterworth filter
    fly.head.pos        = filtfilt(b,a,FlyState{:,2});
  	fly.wing.pos        = filtfilt(b,a,filtfilt(b,a,FlyState{:,3}) - filtfilt(b,a,FlyState{:,4}));
	fly.head.vel        = [diff(fly.head.pos)./fly.Ts ; 0];
	fly.wing.vel        = [diff(fly.wing.pos)./fly.Ts ; 0];
    pat.time            = AI{:,1};
    pat.Ts              = mean(diff(pat.time));
    pat.Fs              = 1/pat.Ts;
 	pat.pos            = AI{:,2};
    
	
 	fly.head.pos        = interp1(fly.time, fly.head.pos , vid.time, 'nearest'); % interpolate pattern y-pos to match fly
 	fly.head.vel        = interp1(fly.time, fly.head.vel , vid.time, 'nearest'); % interpolate pattern y-pos to match fly
 	pat.pos         	= interp1(pat.time, pat.pos      , vid.time, 'nearest'); % interpolate pattern x-pos to match fly
 	
    fly.time            = fly.time      (span);
    fly.head.pos        = fly.head.pos  (span);
	fly.head.vel        = fly.head.vel  (span);
 	pat.pos             = pat.pos       (span);

    TIME                {Index{kk,1}}{Index{kk,3},1}(:,end+1) = fly.time;
    HEAD.Pos            {Index{kk,1}}{Index{kk,3},1}(:,end+1) = fly.head.pos;
  	PAT.Pos            {Index{kk,1}}{Index{kk,3},1}(:,end+1) = pat.pos ;
    
    
    
end






%% Plot head position %%
%---------------------------------------------------------------------------------------------------------------------------------
if showplot.Time
figure (1) ; clf ; hold on
set(gcf,'Color','w')
set(gcf,'Name','Head Position')
set(gcf,'Position',[0 0 Number{1,4}*400 Number{1,3}*400])
movegui(gcf,'center')
for kk = 1:Number{1,1}
    pp = 1;
        for ii = 1:Number{1,3}
            subplot(Number{1,3},1,pp) ; hold on
                xlim([0 20])
                ylim(20*[-1 1])
                if any(pp==(1:Number{1,3}))
                    title(['Amplitude = ' num2str(Unique{1,3}{1}(jj))],'FontSize',12,'fontweight','bold')
                    xlabel('Time (s)','FontSize',12,'fontweight','bold')
                else
                    yticks(0)
                    yticklabels('')  
                    xticks(0)
                    xticklabels('')
                end

                time = TIME{kk}{ii};
                pos  = rad2deg(HEAD.Pos{kk}{ii});
                plot(time,pos)

            pp = pp + 1;
        end
end


