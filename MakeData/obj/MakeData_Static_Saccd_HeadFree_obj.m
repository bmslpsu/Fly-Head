function [] = MakeData_Static_Saccd_HeadFree_obj(Fc)
%% MakeData_Ramp_HeadFree_obj: Reads in all raw trials, transforms data, and saves in organized structure for use with figure functions
%   INPUTS:
%       root    : root directory
%   OUTPUTS:
%       -
%---------------------------------------------------------------------------------------------------------------------------------
% Fc = 30;
filename = ['Ramp_HeadFree_SACCD_filt=' num2str(Fc)];
rootdir = 'H:\EXPERIMENTS\Experiment_Static_SpatFreq';

%% Setup Directories %%
%---------------------------------------------------------------------------------------------------------------------------------
root.daq = rootdir;
root.ang = fullfile(root.daq,'\Vid\HeadAngles\');

% Select files
[D,I,N,U,T,FILES,PATH.ang] = GetFileData(root.ang,false);

PATH.daq = root.daq;

%% Get Data %%
%---------------------------------------------------------------------------------------------------------------------------------
disp('Loading...')
TRIAL = cell(N{1,1},N{1,3});

clear SACD
SACD.Head = [];
SACD.Wing = [];

tt = (0:(1/200):10)';
for kk = 1:N.file
    disp(kk)
    % Load HEAD & DAQ data
	load(fullfile(PATH.daq, FILES{kk}),'data','t_p'); % load pattern x-position
    load(fullfile(PATH.ang, FILES{kk}),'hAngles','t_v'); % load head angles % time arrays
   	%-----------------------------------------------------------------------------------------------------------------------------
    % Get head data
    head.Time = t_v;
    head.Pos = hAngles;
    head.Fc = Fc;
    Head = Fly(head.Pos,head.Time,head.Fc,[],tt); % head object
  	%-----------------------------------------------------------------------------------------------------------------------------
    % Get wing data from DAQ
	wing.f          = medfilt1(100*data(:,6),3); % wing beat frequency [Hz]
    wing.Time       = t_p; % wing time [s]
    wing.Fs         = 1/mean(diff(wing.Time)); % sampling frequency [Hz]
    wing.Fc         = Fc; % cutoff frequency [Hz]
    [b,a]           = butter(2,wing.Fc/(wing.Fs/2)); % butterworth filter
	wing.Left       = filtfilt(b,a,(data(:,4))); % left wing [V]
    wing.Right      = filtfilt(b,a,(data(:,5))); % right wing [V]
    wing.Pos        = wing.Left - wing.Right; % dWBA (L-R) [V]
  	wing.Pos        = wing.Pos - mean(wing.Pos); % subtract mean [V]
   	Wing            = Fly(wing.Pos,t_p,15,[],Head.Time); % wing object
    
    wing.f          = interp1(wing.Time,wing.f,Head.Time);
   	wing.Left     	= interp1(wing.Time,wing.Left,Head.Time);
   	wing.Right     	= interp1(wing.Time,wing.Right,Head.Time);

    Wing.WBF        = wing.f;
    Wing.WBA        = [wing.Left,wing.Right,wing.Left + wing.Right];
    %-----------------------------------------------------------------------------------------------------------------------------
  	% Check WBF & WBA
    if min(wing.f)<150 || mean(wing.f)<200 % check WBF, if too low then don't use trial
        fprintf('Low WBF: Fly %i Trial %i \n',D{kk,1},D{kk,2})
        continue
    elseif any(wing.Left>10.6) || any(wing.Right>10.6)
        fprintf('WBA out of range: Fly %i Trial %i \n',D{kk,1},D{kk,2})
        % continue
    end
    %-----------------------------------------------------------------------------------------------------------------------------
    % Get Saccade Stats   
    [head.SACD,head.thresh,head.count,head.rate,head.SACDRmv] = SacdDetect(Head.X(:,1),Head.Time,2.5,false);
    [wing.SACD,wing.thresh,wing.count,wing.rate,wing.SACDRmv] = SacdDetect(Wing.X(:,1),Wing.Time,1.75,false);
    
    HeadRmv = Fly(head.SACDRmv,Head.Time,[],[],tt);
    WingRmv = Fly(wing.SACDRmv,Wing.Time,[],[],tt);
 	
    head.match = table(head.SACD.Direction*nan);
    mIdx = 1:head.count;
    head.match.Properties.VariableNames = {'Match'};
  	wing.match = table(wing.SACD.Direction*nan);
    wing.match.Properties.VariableNames = {'Match'};
    
    if isnan(head.count)
        head.Rate = table(0);
        wing.Rate = table(0);
    else
        head.Rate = table(nan(head.count,1));
        head.Rate{1,1} = head.rate;
    	wing.Rate = table(nan(wing.count,1));
        wing.Rate{1,1} = wing.rate;
    end
    head.Rate.Properties.VariableNames = {'Rate'};
  	wing.Rate.Properties.VariableNames = {'Rate'};
    
	head.SACD = [head.SACD , head.match];
    wing.SACD = [wing.SACD , wing.match];
	
    I_table = [I(kk,1:3) , rowfun(@(x) abs(x), D(kk,3))];
    I_table.Properties.VariableNames{3} = 'WaveIdx';
    I_table.Properties.VariableNames{4} = 'Wave';
    
    if isnan(head.count)
        head.I_table = I_table;
    else
        head.I_table = repmat(I_table,head.count,1);
    end
    head.SACD = [head.SACD , head.Rate];
    
    if isnan(wing.count)
        wing.I_table = I_table;
    else
        wing.I_table = repmat(I_table,wing.count,1);
    end
    
	head.SACD = [head.I_table , head.SACD];
  	wing.SACD = [wing.I_table , wing.SACD];
    
    SACD.Head = [SACD.Head ; head.SACD];
    SACD.Wing = [SACD.Wing ; wing.SACD];
        
 	%-----------------------------------------------------------------------------------------------------------------------------
    % Store objects in cells
    vars = {[],Head,Wing,HeadRmv,WingRmv};
	qq = size(TRIAL{I{kk,1},I{kk,3}},1);
    for ww = 1:length(vars)
        TRIAL{I{kk,1},I{kk,3}}{qq+1,ww} = vars{ww};
    end

%     pause()
%     close all
end

clear jj ii kk pp qq ww n a b  t_v hAngles data head wing pat tt I_table Dir loop Saccade Interval Stimulus Error IntError...
    Head Pat Wing  vars root t_p var1 var2 Err_table pos_err vel_err pos_int_err vel_int_err stim_pos matchFlag emptyFlag
                        
%% Fly Statistics %%
%---------------------------------------------------------------------------------------------------------------------------------
FLY = cell(N{1,3},1);
for kk = 1:N{1,1}
    for jj = 1:N{1,3}
        for ii = 1:size(TRIAL{kk,jj},2)
            FLY{jj}{kk,ii} = FlyStats(TRIAL{kk,jj}(:,ii));
        end
    end
end
clear kk jj ii
%% Grand Statistics %%
%---------------------------------------------------------------------------------------------------------------------------------
GRAND = cell(size(FLY,2),size(FLY{1},2));
for jj = 1:N{1,3}
    for ii = 1:size(FLY{jj},2)
        GRAND{jj,ii} = GrandStats(FLY{jj}(:,ii));
    end
end
clear jj ii
%% SAVE %%
%---------------------------------------------------------------------------------------------------------------------------------
disp('Saving...')
save(['H:\DATA\Rigid_Data\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
    'SACD','TRIAL','FLY','GRAND','D','I','U','N','T','-v7.3')
disp('SAVING DONE')
end