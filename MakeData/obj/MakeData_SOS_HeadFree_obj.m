function [] = MakeData_SOS_HeadFree_obj(rootdir)
%% MakeData_SOS_HeadFree_obj: Reads in all raw trials, transforms data, and saves in organized structure for use with figure functions
%   INPUTS:
%       rootdir    	:   root directory
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
rootdir = 'H:\EXPERIMENTS\Experiment_SOS_v2';
filename = 'SOS_HeadFree_DATA';
%---------------------------------------------------------------------------------------------------------------------------------
%% Setup Directories %%
%---------------------------------------------------------------------------------------------------------------------------------
root.daq = rootdir;
% root.ang = fullfile(root.daq,'\Vid\Angles\');
root.ang = fullfile(root.daq,'\Angles\');

% Select files
[FILES, PATH.ang] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root.ang, 'MultiSelect','on');
FILES = cellstr(FILES)';

PATH.daq = root.daq;

[D,I,N,U,T] = GetFileData(FILES);

clear rootdir
%% Get Data %%
%---------------------------------------------------------------------------------------------------------------------------------
% IOFreq = 0.1*round(linspace(0.1,8,10)/0.1)';
IOFreq = [1, 3.1, 5.3, 7.4, 9.6];
disp('Loading...')
ALL 	= cell([N{1,end},10]); % cell array to store all data objects
TRIAL  	= cell(N{1,1},1);
n.catg  = size(N,2) - 1;
pp = 0;
span = 1:2000;
tt = linspace(0,20,2000)';
for kk = 1:N{1,end}
    disp(kk)
    % Load HEAD & DAQ data
    data = [];
	load(fullfile(PATH.daq, FILES{kk}),'data','t_p'); % load pattern x-position
    load(fullfile(PATH.ang, FILES{kk}),'hAngles','t_v'); % load head angles % time arrays
	%-----------------------------------------------------------------------------------------------------------------------------
    % Get head data
    head.Time = t_v(span);
    head.Pos = hAngles(span) - mean(hAngles(span));
    Head = Fly(head.Pos,head.Time,40,IOFreq); % head object
	%-----------------------------------------------------------------------------------------------------------------------------
    % Get wing data from DAQ
	wing.f          = medfilt1(100*data(:,6),3); % wing beat frequency [Hz]
    wing.Time       = t_p; % wing time [s]
    wing.Fs         = 1/mean(diff(wing.Time)); % sampling frequency [Hz]
    wing.Fc         = 20; % cutoff frequency [Hz]
    [b,a]           = butter(2,wing.Fc/(wing.Fs/2)); % butterworth filter
	wing.Left       = filtfilt(b,a,(data(:,4))); % left wing [V]
    wing.Right      = filtfilt(b,a,(data(:,5))); % right wing [V]
    wing.Pos        = wing.Left - wing.Right; % dWBA (L-R) [V]
  	wing.Pos        = wing.Pos - mean(wing.Pos); % subtract mean [V]
   	Wing            = Fly(wing.Pos,t_p,40,IOFreq,Head.Time); % wing object
    
    wing.f          = interp1(wing.Time,wing.f,Head.Time);
   	wing.Left     	= interp1(wing.Time,wing.Left,Head.Time);
   	wing.Right     	= interp1(wing.Time,wing.Right,Head.Time);

    Wing.WBF        = wing.f;
    Wing.WBA        = [wing.Left,wing.Right,wing.Left + wing.Right];
    %-----------------------------------------------------------------------------------------------------------------------------
    % Check WBF & WBA
    if min(wing.f)<150 || mean(wing.f)<180 % check WBF, if too low then don't use trial
        fprintf('Low WBF: Fly %i Trial %i \n',D{kk,1},D{kk,2})
        continue
    elseif any(wing.Left>10.6) || any(wing.Right>10.6)
        fprintf('WBA out of range: Fly %i Trial %i \n',D{kk,1},D{kk,2})
        continue
    else
        pp = pp + 1; % set next index to store data
    end
	%-----------------------------------------------------------------------------------------------------------------------------
	% Get pattern data from DAQ
    pat.Time	= t_p;
    pat.Pos 	= panel2deg(data(:,2));  % pattern x-pos: subtract mean and convert to deg [deg]  
    pat.Pos  	= FitPanel(pat.Pos,pat.Time,tt); % fit panel data
 	Pat      	= Fly(pat.Pos,Head.Time,0.4*Head.Fs,IOFreq); % pattern object
	%-----------------------------------------------------------------------------------------------------------------------------
 	% Calculate error between head & pattern
    head.Err = Pat.X(:,1) - Head.X(:,1); % calculate position error between head & pattern [deg]
    Err = Fly(head.Err,Head.Time,0.4*Head.Fs,IOFreq); % error object
    %-----------------------------------------------------------------------------------------------------------------------------
    % Calculate iput-output relationships
    pat2head    = IO_Class(Pat,Head);
    err2wing    = IO_Class(Err,Wing);
	head2wing   = IO_Class(Head,Wing);
    pat2wing    = IO_Class(Pat,Wing);
    %-----------------------------------------------------------------------------------------------------------------------------
    % Store objects in cells
    for jj = 1:n.catg
        ALL{pp,jj} = I{kk,jj};
    end

	vars = {Pat,Head,Wing,Err,pat2head,err2wing,head2wing,pat2wing};
    for jj = 1:length(vars)
        ALL{pp,n.catg+jj} = vars{jj};
    end
    
	qq = size(TRIAL{I{kk,1}},1);
    for ww = 1:length(vars)
        TRIAL{I{kk,1}}{qq+1,ww} = vars{ww};
    end
end

ALL( all(cellfun(@isempty, ALL),2), : ) = []; % get rid of emtpty rows becuase of low WBF

clear jj ii kk pp qq ww n a b spant_p t_v hAngles data head wing pat bode tt ...
    Head Pat Wing Err pat2head pat2wing err2wing head2wing vars root t_p span
disp('LOADING DONE')

%% Fly Statistics %%
%---------------------------------------------------------------------------------------------------------------------------------
FLY = cell(N{1,1},size(TRIAL{1},2));
for kk = 1:N{1,1}
    for ii = 1:size(TRIAL{kk},2)
        FLY{kk,ii} = FlyStats(TRIAL{kk}(:,ii));
    end
end
clear kk ii
%% Grand Statistics %%
%---------------------------------------------------------------------------------------------------------------------------------
GRAND = cell(1,size(FLY,2));
for ii = 1:size(FLY,2)
    GRAND{ii} = GrandStats(FLY(:,ii));
end
clear ii

%% Replay Stimulus %%
%---------------------------------------------------------------------------------------------------------------------------------
errIdx = 4;
xIdx = 1;

Time  	= GRAND{1}.Mean{1}{5}(:,1);
RefMean	= GRAND{1}.Mean{1}{6}(:,xIdx);
ErrMean = GRAND{errIdx}.Mean{1}{6}(:,xIdx);
ErrSTD	= GRAND{errIdx}.STD{1}{6}(:,xIdx);

Error = [];
dev = [];
pp = 1;
for kk = 1:N.fly
    for ii = 1:size(TRIAL{kk},1)
        Error(:,pp) = TRIAL{kk}{ii,errIdx}.X(:,xIdx);
              
        pp = pp + 1;
    end
end

dev = abs(ErrMean - Error);
devMean = mean(dev,1);
[~,minIdx] = min(devMean);

TEST = Fly(Error(:,minIdx),Time,[],IOFreq);

FIG = figure (1) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 8 4];
movegui(FIG, 'center')

subplot(2,1,1) ; hold on
xlabel('Time')
ylabel('Angle (deg)')
% plot(Time,RefMean,'b')
[~,h.ErrMean] = PlotPatch(ErrMean, ErrSTD, Time, 3, 1, 'k' ,[0.4 0.4 0.6], 0.4 , 2);
h.ErrTrial = plot(Time,Error(:,minIdx),'r');

subplot(2,1,2) ; hold on
xlim([0 10])
xlabel('Frequency (Hz)')
ylabel('Mag (deg)')
plot(TRIAL{1}{1,1}.Fv, TRIAL{1}{1,1}.Mag(:,xIdx), 'b', 'LineWidth', 1)
plot(TEST.Fv, TEST.Mag(:,1), 'r', 'LineWidth', 1)

plot(TRIAL{1}{1,1}.IOFreq, TRIAL{1}{1,1}.IOMag(:,xIdx), '-b*', 'LineWidth', 1)
plot(TEST.IOFreq, TEST.IOMag(:,1), '-r*', 'LineWidth', 1)

legend('Reference','Error')

func = 3.75*deg2panel(Error(:,minIdx));
Fs = 200;
TT = round(Time(end));
tt = 0:(1/Fs):TT;
func = interp1(Time, func , tt, 'linear','extrap'); % interpolate to match new time

subplot(2,1,1) ; hold on
h.func = plot(tt, func, 'g', 'LineWidth', 1);
legend([h.ErrMean,h.ErrTrial,h.func], 'Mean Error','Selected Trial Error','Function')

% Name file
strFreq = '';
for kk = 1:length(IOFreq)
   strFreq = [strFreq  num2str(IOFreq(kk)) '_'];
end
strFreq(end) = [];

fname = sprintf(['position_function_SOS_REPLAY_fs_%1.1f_T_%1.1f_freq_' strFreq '.mat'],Fs,TT);
targetDir = 'C:\Users\boc5244\Documents\GitHub\Arena\Functions';
save(fullfile(targetDir,fname), 'func');

%% SAVE %%
%---------------------------------------------------------------------------------------------------------------------------------
disp('Saving...')
save(['H:\DATA\Rigid_Data\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
    'TRIAL','FLY','GRAND','D','I','U','N','T','-v7.3')
disp('SAVING DONE')
end