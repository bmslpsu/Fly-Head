function [] = MakeData_Sine_HeadFree_Intervals(rootdir,Amp)
%% MakeData_Sine_HeadFree: Reads in all raw trials, transforms data, and saves in organized structure for use with figure functions
%   INPUTS:
%       root  	: root directory
%       Amp   	: sinusoid amplitude, set's directory
%   OUTPUTS:
%       PAT     : pattern structure
%       WING   	: wings structure
%       HEAD  	: head structure
%---------------------------------------------------------------------------------------------------------------------------------
% EXAMPLE INPUT %
Amp = 15;
rootdir = 'F:\EXPERIMENTS\Experiment_Sinusoid\';
%---------------------------------------------------------------------------------------------------------------------------------
%% Setup Directories %%
%---------------------------------------------------------------------------------------------------------------------------------
% filename = ['Sine_HeadFree_Interval' num2str(Amp) '_DATA_' ] ;
root.daq = [rootdir num2str(Amp) '\'];
root.ang = [root.daq '\Vid\Angles\'];

% Select files
[FILES, PATH.ang] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root.ang, 'MultiSelect','on');
FILES = FILES';

PATH.daq = root.daq;

[D,I,N,U] = GetFileData(FILES,'Fly','Trial','Freq');

clear rootdir
%% Get Data %%
%---------------------------------------------------------------------------------------------------------------------------------
disp('Loading...')
% Preallocate data cells
PAT             = cell(N{1,1},1);
HEAD            = cell(N{1,1},1);
WING            = cell(N{1,1},1);
ERR             = cell(N{1,1},1);
IO.pat2head     = cell(N{1,1},1);
IO.err2wing     = cell(N{1,1},1);
IO.head2wing    = cell(N{1,1},1);
for kk = 1:N{1,1}
    PAT{kk}             = cell(N{1,3},1);
    HEAD{kk}            = cell(N{1,3},1);
    WING{kk}            = cell(N{1,3},1);
    ERR{kk}             = cell(N{1,3},1);
    IO.pat2head{kk}     = cell(N{1,3},1);
    IO.err2wing{kk}     = cell(N{1,3},1);
    IO.head2wing{kk}    = cell(N{1,3},1);
end
nInt = 5;
% Store data in cells
count = 0;
pp = 0;
for kk = 1:N{1,4}
    disp(kk)
    % Load HEAD & DAQ data
    data = [];
	load([PATH.daq   FILES{kk}],'data','t_p'); % load pattern x-position
    load([PATH.ang   FILES{kk}],'hAngles','t_v'); % load head angles % time arrays
    Freq = U{1,3}{1}(I{:,3}(kk)); % frequency input for trial
    % Check WBF
	wing.f = 100*(data(:,6)); % wing beat frequency
    if min(wing.f)<150 || mean(wing.f)<180 % check WBF, if too low then don't use trial
        fprintf('Low WBF: Fly %i Trial %i \n',FD.Fly(kk),FD.Trial(kk))
        continue
    else
        pp = pp + 1; % set next index to store data
    end
    % Get head data
    head.Time = t_v;
    head.Pos = hAngles - mean(hAngles);
    % Get wing data from DAQ
    wing.Time       = t_p; % wing time [s]
    wing.Fs         = 1/mean(diff(wing.Time)); % sampling frequency [Hz]
    wing.Fc         = 20; % cutoff frequency [Hz]
    [b,a]           = butter(2,wing.Fc/(wing.Fs/2)); % butterworth filter
    wing.Left       = filtfilt(b,a,(data(:,4))); % left wing [V]
    wing.Right      = filtfilt(b,a,(data(:,5))); % right wing [V]
    wing.Pos        = wing.Left - wing.Right; % dWBA (L-R) [V]
    wing.Pos        = wing.Pos - mean(wing.Pos); % subtract mean [V]
    
    wing.Pos = interp1(wing.Time, wing.Pos , head.Time, 'nearest'); % interpolate to match new time
    
    % Get pattern data from DAQ
    pat.Time	= t_p;
    pat.Pos 	= panel2deg(data(:,2));  % pattern x-pos: subtract mean and convert to deg [deg]  
    pat.Pos  	= FitPanel(pat.Pos,pat.Time,head.Time,false,false); % fit panel data
    % Calculate error between head & pattern
    %-----------------------------------------------------------------------------------------------------------------------------
    qq = 0;
    interval = round(length(head.Time)/nInt);
 	idx = size(PAT{I{:,1}(kk)}{I{:,3}(kk)},1);
    for jj = 1:nInt
        span = (qq+1):qq+interval;
        % Objects
        Head    = Fly(head.Pos(span),head.Time(span),20); % head object
        Head    = IO_Freq(Head,Freq);
        
        Wing    = Fly(wing.Pos(span),pat.Time(span),20); % wing object
        Wing    = IO_Freq(Wing,Freq);
        
        Pat     = Fly(pat.Pos(span),head.Time(span),0.4*Head.Fs); % pattern object
        Pat     = IO_Freq(Pat,Freq);
        
        head.Err = Pat.X(:,1) - Head.X(:,1); % calculate position error between head & pattern [deg]

        Err     = Fly(head.Err,head.Time(span),0.4*Head.Fs); % error object
        Err     = IO_Freq(Err,Freq);

        % Calculate iput-output relationships
        pat2head    = IO_Class(Pat,Head);
        err2wing    = IO_Class(Err,Head);
        head2wing   = IO_Class(Head,Wing);
        
        % CELLS
        
        PAT             {I{:,1}(kk)}{I{:,3}(kk)}{idx+1,jj} = Pat;
        HEAD            {I{:,1}(kk)}{I{:,3}(kk)}{idx+1,jj} = Head;
        WING            {I{:,1}(kk)}{I{:,3}(kk)}{idx+1,jj} = Wing;
        ERR             {I{:,1}(kk)}{I{:,3}(kk)}{idx+1,jj} = Err;
        IO.pat2head  	{I{:,1}(kk)}{I{:,3}(kk)}{idx+1,jj} = pat2head;
        IO.err2wing 	{I{:,1}(kk)}{I{:,3}(kk)}{idx+1,jj} = err2wing;
        IO.head2wing  	{I{:,1}(kk)}{I{:,3}(kk)}{idx+1,jj} = head2wing;
        
        qq = qq + length(span);
    end
    %-----------------------------------------------------------------------------------------------------------------------------
end
disp('LOADING DONE')
disp('Bad Trial:')
disp(count)
pause(3)
clear jj kk pp qq a b t_p t_v hAngles data head wing pat bode tt root pat2head err2wing head2wing span interval idx ...
    Head Wing Pat Err count Amp
%%
%---------------------------------------------------------------------------------------------------------------------------------
FIG = figure (1) ; clf ; hold on
FIG.Color = 'w';
xlabel('Frequency (Hz)')
ylabel('Phase (deg)')
ALL = cell(N{1,3},nInt);
cFlow = jet(nInt);
for kk = 1:N{1,1}
    for jj = 1:N{1,3}
        for ii = 1:size(IO.pat2head{kk}{jj},1)
            for ww = 1:size(IO.pat2head{kk}{jj},2)
                ALL{jj,ww}(end+1,1) = IO.pat2head{kk}{jj}{ii,ww}.IOBodePhaseDiff(:,1);               
                plot(IO.pat2head{kk}{jj}{ii,ww}.IOFreq(:,1),rad2deg(IO.pat2head{kk}{jj}{ii,ww}.IOBodePhaseDiff(:,1)),...
                   '-+','Color',cFlow(ww,:))
                xlim([0 4])
                ylim([-100 300])
            end
        end
    end
end

ALL_Mean = cell2mat(cellfun(@(x) circ_mean(x,[],1),ALL,'UniformOutput',false));

% figure (1) ; clf ; hold on
% xlim([0 13])
for ww = 1:nInt
    plot(U{1,3}{1},rad2deg(ALL_Mean(:,ww)),'Color',cFlow(ww,:),'LineWidth',3)
end

%%
%---------------------------------------------------------------------------------------------------------------------------------
cFlow = jet(N{1,3});
figure (2) ; clf ; hold on
xlabel('Interval')
ylabel(['Phase (' char(176) ')'])
% xlim([1 nInt])
xticks(1:nInt)
for ww = 1:N{1,3}
    plot(1:nInt,rad2deg(ALL_Mean(ww,:)),'Color',cFlow(ww,:),'LineWidth',3)
end

legList = cell(N{1,3},1);
for kk = 1:N{1,3}
   legList{kk} = [num2str(U{1,3}{1}(kk)) ' Hz'];    
end
legend(legList)

% for kk = 1:N{1,1}
%     for jj = 1:N{1,3}
%         for ii = 1:size(IO.pat2head{kk}{jj},1)
%             for ww = 1:size(IO.pat2head{kk}{jj},2)
%                 plot(ww,rad2deg(IO.pat2head{kk}{jj}{ii,ww}.IOBodePhaseDiff(:,1)),...
%                    '-+','Color',cFlow(jj,:))
%             end
%         end
%     end
% end




%% Save ouputs as structure %%
%---------------------------------------------------------------------------------------------------------------------------------
% disp('Saving...')
% save([rootdir 'DATA\' filename '.mat'],'PAT','WING','HEAD','BODE','CROSS','D','I','N','U')
% disp('DONE')
end