

clear
showplot.Time = 0;
showplot.Sacd = 0;

%% Setup Directories %%
%---------------------------------------------------------------------------------------------------------------------------------
root.pat = 'H:\EXPERIMENTS\Experiment_Ramp\';
root.head = [root.pat '\Vid\Angles\'];

% Select files
[FILES, PATH.head] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root.head, 'MultiSelect','on');
FILES = FILES';

PATH.pat = uigetdir(root.pat);
PATH.pat = [PATH.pat '\'];

%% Process File Data %%
%---------------------------------------------------------------------------------------------------------------------------------
% clearvars -except FILES PATH Amp
% close all
clc
% Read in data from head file names: [Fly=FLy#, Trial=trial# for each fly]nTrial = length(FILES);     % total # of trials
nTrial = length(FILES); % total # of trials
Fly = zeros(nTrial,1); Trial = zeros(nTrial,1); Vel = zeros(nTrial,1); HEAD.FileCells = cell(nTrial,6);% preallocate arrays
for jj = 1:nTrial
    temp = textscan(char(FILES{jj}), '%s', 'delimiter', '_.'); temp = temp{1} ; % read individual strings into temp variable
    HEAD.FileCells(jj,:) = {temp{1} temp{2} temp{3} temp{4} temp{5} temp{6}}; % separate strings into six rows with one item in each cell
    Fly(jj,1) = str2double(temp{2}); % store fly #
    Trial(jj,1) = str2double(temp{4}); % store trial #
	Vel(jj,1) = str2double(temp{6});   
end
clear temp
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

% Make indexing array for velocity
uVel = sort(unique(Vel)); % find all unique velocities
nVel = length(uVel); % # of unique velocities
idxVel = zeros(nTrial,1);
for kk = 1:nVel
   idx = find(Vel == uVel(kk));
   idxVel(idx) = kk; % velocity index
end
uVel = 3.75*uVel;
fprintf('Total Flies''s: %i \n',nFly) ; fprintf('Total Trials''s: %i \n',nTrial)
T = cell2table(trialFly,'VariableNames',{'Fly','Trials'});
disp(T)

%% Get Data %%
%---------------------------------------------------------------------------------------------------------------------------------
disp('Loading Data...')
clear HEAD PAT WING head pat wings Amplitude
% Preallocate data cells %
for kk = 1:nFly % # of flys
    HEAD.Time{kk,1}         = cell(nVel,1);
    HEAD.Pos{kk,1}          = cell(nVel,1);
    HEAD.Vel{kk,1}          = cell(nVel,1);
    HEAD.VelMed{kk,1}       = cell(nVel,1);
    HEAD.VelSTD{kk,1}       = cell(nVel,1);
    HEAD.Err.Pos{kk,1}      = cell(nVel,1);
    HEAD.Err.Vel{kk,1}      = cell(nVel,1);
    HEAD.ErrSum.Pos{kk,1} 	= cell(nVel,1);
    HEAD.ErrSum.Vel{kk,1} 	= cell(nVel,1);
    HEAD.Freq{kk,1}         = cell(nVel,1);
	HEAD.Mag{kk,1}          = cell(nVel,1);

    PAT.Time{kk,1}          = cell(nVel,1);
    PAT.Pos{kk,1}           = cell(nVel,1);
    PAT.Vel{kk,1}           = cell(nVel,1);
end
% Save data in cells %
for kk = 1:nTrial
    clear head pat data t_v t_p hAngles % clear temporary variables
    %-----------------------------------------------------------------------------------------------------------------------------
    % Load head & DAQ data %
    load([PATH.head  FILES{kk}]); % load head angles % time arrays
	load([PATH.pat   FILES{kk}]); % load pattern x-position
	%-----------------------------------------------------------------------------------------------------------------------------
	% Check WBF %
    wing.f = 100*(data(:,6)); % wing beat frequency
    if min(wing.f)<180
        fprintf('Low WBF: Fly %i Trial %i \n',Fly(kk),Trial(kk))
    end
	%-----------------------------------------------------------------------------------------------------------------------------
	% Get head data %
	head.Time = t_v; % store head time vector
    
	% Setup filter for head %
    head.Fs = 1/mean(diff(head.Time)); % sampling frequency [Hz]
    head.Fc = 20; % cutoff frequency [Hz]
    [b,a] = butter(2,head.Fc/(head.Fs/2)); % butterworth filter
    
    head.Pos = filtfilt(b,a,hAngles); % filter head position [deg]
    head.Pos = head.Pos; % subtract DC component
    head.Vel = filtfilt(b,a,diff(head.Pos)./(1/head.Fs)); % filtered velocity [deg/s]
    head.Vel = [head.Vel ; head.Vel(end)]; % filtered velocity [deg/s]
    head.VelMed = median(abs(head.Vel)); % velocity median
    head.VelSTD = std(abs(head.Vel)); % velcoity STD
    
    [head.Freq, head.Mag , ~] = FFT(head.Time,head.Pos);
    %-----------------------------------------------------------------------------------------------------------------------------
    % Get pattern data from DAQ %
%     pat.Time = t_p; % pattern time
    pat.Time = linspace(0,10,length(head.Time))';
    head.Time = pat.Time;
%     pat.Pos = data(:,2); % pattern position [deg]
    pat.Fs = 1/mean(diff(pat.Time)); % pattern sampling frequency [Hz]
    m = 3.75*Vel(kk);
    pat.Pos = pat.Time*m; 
    pat.Vel = 3.75*Vel(kk)*ones(length(head.Time),1);
 	%-----------------------------------------------------------------------------------------------------------------------------
    % Calculate Head Error %
    head.Err.Pos = head.Pos - pat.Pos; % position error
    head.Err.Vel = head.Vel - pat.Vel; % velcoity error
  	head.ErrSum.Pos = trapz(head.Time,head.Err.Pos); % position error sum
  	head.ErrSum.Vel = trapz(head.Time,head.Err.Vel); % position velcoity sum
	%-----------------------------------------------------------------------------------------------------------------------------
    % Store data in cells %
    % Head
	HEAD.Time       {idxFly(kk),1}{idxVel(kk),1}(:,end+1) = head.Time;      % head time
	HEAD.Pos        {idxFly(kk),1}{idxVel(kk),1}(:,end+1) = head.Pos;       % head position
    HEAD.Vel        {idxFly(kk),1}{idxVel(kk),1}(:,end+1) = head.Vel;       % head velocity
    HEAD.VelMed     {idxFly(kk),1}{idxVel(kk),1}(:,end+1) = head.VelMed;    % mean head velocity
    HEAD.VelSTD     {idxFly(kk),1}{idxVel(kk),1}(:,end+1) = head.VelSTD;    % STD head velocity
    HEAD.Err.Pos    {idxFly(kk),1}{idxVel(kk),1}(:,end+1) = head.Err.Pos;   % head position error
    HEAD.Err.Vel    {idxFly(kk),1}{idxVel(kk),1}(:,end+1) = head.Err.Vel;   % head velocity error
    HEAD.ErrSum.Pos	{idxFly(kk),1}{idxVel(kk),1}(:,end+1) = head.Err.Pos;   % head position error sum
    HEAD.ErrSum.Vel	{idxFly(kk),1}{idxVel(kk),1}(:,end+1) = head.Err.Vel;   % head velocity error sum
    HEAD.Freq       {idxFly(kk),1}{idxVel(kk),1}(:,end+1) = head.Freq;    	% head frequency vector
    HEAD.Mag        {idxFly(kk),1}{idxVel(kk),1}(:,end+1) = head.Mag;       % head position FFT

    % Pattern
	PAT.Time  	{idxFly(kk),1}{idxVel(kk),1}(:,end+1) = pat.Time;   % pattern time
	PAT.Pos    	{idxFly(kk),1}{idxVel(kk),1}(:,end+1) = pat.Pos;    % pattern position
	PAT.Vel    	{idxFly(kk),1}{idxVel(kk),1}(:,end+1) = pat.Vel;    % pattern velocity
    %-----------------------------------------------------------------------------------------------------------------------------
    if showplot.Time
    rows = 4;
    figure (100)
        subplot(ceil(nTrial/rows),rows,kk) ; hold on
            title(['Fly ' num2str(Fly(kk)) ' Vel ' num2str(Vel(kk))])
            plot(pat.Time,pat.Pos,'k')
            plot(head.Time,head.Pos,'b')
%             plot(head.Time,head.Err.Pos,'r')
            box on
            hold off
            xlim([0 10])
%             ylim([-1000 1000])
   
    figure (101)
        subplot(ceil(nTrial/rows),rows,kk) ; hold on
            title(['Fly ' num2str(Fly(kk)) ' Vel ' num2str(Vel(kk))])
            plot(head.Time,pat.Vel,'k','LineWidth',3)
            plot(head.Time,head.Vel,'b')
            plot(head.Time,head.Err.Vel,'-r')
            box on
            hold off
            xlim([0 10])
%             ylim([-1000 1000])
    end
    %-----------------------------------------------------------------------------------------------------------------------------
end
disp('DONE')

%% Example Trial %%
%---------------------------------------------------------------------------------------------------------------------------------
figure (22) ; clf
subplot(2,1,1) ; hold on ; box on
set(gca,'FontSize',28)
plot(pat.Time,pat.Pos,'k','LineWidth',5)
plot(head.Time,head.Pos,'b','LineWidth',5)
xlim([0 3])
xticks([])
yticks([0 50 100])
ylabel(['Angle (' char(176) ')'])
legend('Visual Scene','Head')

subplot(2,1,2) ; hold on ; box on
set(gca,'FontSize',28)
plot(pat.Time,pat.Vel,'k','LineWidth',5)
plot(head.Time,head.Vel,'b','LineWidth',5)
xlim([0 3])
xticks([0 1 2 3])
ylim([-400 400])
% yticks([-])
xlabel('Time (s)')
ylabel(['Velcoity (' char(176) '/s)'])

%% FFT of HEAD Position %%
%---------------------------------------------------------------------------------------------------------------------------------
figure (23) ; clf ; hold on
colorIdx = {'b','g','r','r','g','b'};
pp = 1;
for jj = [4,3,5,2,6,1]
	figure (23) ; subplot(3,2,pp) ; hold on ; title(['Vel = ' num2str(uVel(jj))])
    xlabel('Frequency (Hz)') ; ylabel('|P1|') ; box on
    for kk = 1:nFly
        for ii = 1:size(HEAD.Freq{kk}{jj},2)
            plot(HEAD.Freq{kk}{jj}(:,ii),HEAD.Mag{kk}{jj}(:,ii),colorIdx{jj})
            xlim([4 20])
            ylim([0 1])
        end
    end
    pp = pp + 1;
end

%% Saccade Detection %%
%---------------------------------------------------------------------------------------------------------------------------------
disp('Saccade Detection...')
clear Sacs pks 
if showplot.Sacd
    figure (34) ; clf
    figure (35) ; clf
end
for kk = 1:nFly
    HEAD.SACD.INTER.Time        {kk,1} = cell(nVel,1);
    HEAD.SACD.INTER.Pos         {kk,1} = cell(nVel,1);
    HEAD.SACD.INTER.Vel         {kk,1} = cell(nVel,1);
    HEAD.SACD.INTER.PosErr      {kk,1} = cell(nVel,1);
    HEAD.SACD.INTER.VelErr      {kk,1} = cell(nVel,1);
    HEAD.SACD.INTER.IntPosErr 	{kk,1} = cell(nVel,1);
    HEAD.SACD.INTER.IntVelErr 	{kk,1} = cell(nVel,1);
    
    PAT.SACD.INTER.Pos    {kk,1} = cell(nVel,1);
    PAT.SACD.INTER.Vel    {kk,1} = cell(nVel,1);
end

HEAD.SACD.MAIN = double.empty(0,22);

pp = 1;
rows = 5;
for kk = 1:nFly
    for jj = [4 3 5 2 6 1]
        for ii = 1:size(HEAD.Pos{kk}{jj},2)
            SacdThresh = HEAD.VelMed{kk}{jj}(:,ii) + 3*HEAD.VelSTD{kk}{jj}(:,ii); % threshold for saccade detetcion (by trial)
            
            TimeData = HEAD.Time{kk}{jj}(:,ii); % get time data
            VelData = HEAD.Vel{kk}{jj}(:,ii); % get velocity data
            AbsVelData = abs(VelData); % absolute value for velocity
            SacdVel = AbsVelData; SacdVel(SacdVel<SacdThresh) = 0; % set all data below threshold to 0
            PosData = HEAD.Pos{kk}{jj}(:,ii); % get position data
            PatData = PAT.Pos{kk}{jj}(:,ii); % get pattern data
            PatVelData = PAT.Vel{kk}{jj}(:,ii); % get pattern data
         	Fs = 1/mean(diff(TimeData)); % sampling frequency [Hz]
                       
             % Find local maxima
            [pks, pIdx] = findpeaks(SacdVel,'MINPEAKDISTANCE',25);
            step = 30; % window length (in samples)
            ta = 30; % adaptation time (in samples)
            [I] = find(pIdx > (ta) & pIdx < ((length(TimeData) - (step+1)))); % ignore saccades at beginning and end
            pIdx = pIdx(I);
            pks = pks(I);
            clear Sacs
            Sacs = nan(length(pIdx),16);
            mm = 1;
            preWin = 15;
            if ~isempty(pks)
                for ww = 1:length(pIdx)
                    % FLY TRIAL INFO
                    Sacs(mm,1) 	= kk;               % Fly # (normalized)
                    Sacs(mm,2)  = abs(uVel(jj));    % Speed
                    Sacs(mm,3)  = sign(uVel(jj));  	% Direction
                    Sacs(mm,4)  = ii;               % Trial #
                    
                   % DURATION: define interval as 1/4 peak amplitude
                    sIdx(mm) = find(AbsVelData(1:pIdx(ww)) <= pks(ww)/4,1,'last'); % saccade start index               
                    Eind = find(AbsVelData <= pks(ww)/4); % all values below 1/4 peak
                    Es = find(Eind > pIdx(ww),1,'first'); % first value after the start index is the end idex
                    if ~isempty(Es) % make sure data did not start above threshold
                        eIdx(mm) = Eind(Es); % saccade end index
                    end
                    
                 	% TIME
                    Sacs(mm,5) = TimeData(sIdx(mm)); % start time
                    Sacs(mm,6) = TimeData(pIdx(ww)); % peak time
                    Sacs(mm,7) = TimeData(eIdx(mm)); % end time
                    Sacs(mm,8) = 1000*(TimeData(eIdx(mm)) - TimeData(sIdx(mm))); % duration

                    % POSITION & AMPLITUDE
                    Sacs(mm,9)  = PosData(sIdx(mm));  % start position
                    Sacs(mm,10) = PosData(pIdx(ww));  % peak position
                    Sacs(mm,11) = PosData(eIdx(mm));  % end position
                    Sacs(mm,12) = abs(PosData(sIdx(mm)) - PosData(eIdx(mm))); % amplitude

                    % VELOCITY
                    Sacs(mm,13) = AbsVelData(sIdx(mm)); % start velocity
                    Sacs(mm,14) = pks(ww);              % abs peak velocity
                    Sacs(mm,15) = AbsVelData(eIdx(mm)); % end velocity

                    % MAX PRE-SACCADE VELOCITY
                    if ww == 1
                        if pIdx(ww)<=preWin
                            Sacs(mm,16) = NaN;
                        else
                            Sacs(mm,16) = abs(median(VelData(pIdx(ww)-preWin:pIdx(ww))));
                        end
                    elseif (abs((pIdx(ww-1)-pIdx(ww))) > preWin) % make sure fly stops before next saccade
                        Sacs(mm,16) = abs(median(VelData(sIdx(mm)-preWin:sIdx(mm))));
                    else
                        Sacs(mm,16) = NaN;
                    end  

                    % INTER-SACCADE INTERVAL
                    if ww == 1
                        interval = sIdx(mm)-1;
                    else
                        interval = (sIdx(mm) - eIdx(mm-1));
                    end
                    
                    % HEAD INTER-SACCADE: TIME,POSITION,VELOCIY
                    headTime  = TimeData(sIdx(mm)-interval:sIdx(mm));   % head raw head time
                    headPos   = PosData(sIdx(mm)-interval:sIdx(mm));    % head raw position
                    headVel   = VelData(sIdx(mm)-interval:sIdx(mm));    % head raw velocity
                   	NheadTime = headTime - headTime(1);                 % head normalized time
                    NheadPos  = headPos - headPos(1);                   % head normalized position
                    NheadVel  = headVel - headVel(1);                   % head normalized velocity

                    % PAT INTER-SACCADE: POSITION,VELOCIY
                    patPos   = PatData(sIdx(mm)-interval:sIdx(mm));     % pat raw position
                    patVel   = PatVelData(sIdx(mm)-interval:sIdx(mm));  % pat raw velocity
                    NpatPos  = patPos - patPos(1);                      % pat normalized position
                    NpatVel  = patVel - patVel(1);                      % pat normalized velocity
                    
                    % ERROR INTER-SACCADE: POSITION,VELOCIY
                    errPos   = patPos  -  headPos;    % raw position error 
                    errVel   = patVel  -  headVel;    % raw velocity error 
                    NerrPos  = NpatPos -  NheadPos;   % normalized position error
                    NerrVel  = NpatVel -  NheadVel;   % normalized velocity error
                    
                    Sacs(mm,17) = mean(NerrPos(end-5));
                    Sacs(mm,18) = mean(errVel(end-5));
                    
                	% INTEGRATED ERROR INTER-SACCADE: POSITION,VELOCIY
                    INTerrPos  = cumtrapz(NheadTime,NerrPos);	% integrated position error 
                    INTerrVel  = cumtrapz(NheadTime,NerrVel);	% integrated velocity error 
                    
                    Sacs(mm,19) = mean(INTerrPos(end-5));
                    Sacs(mm,20) = mean(INTerrVel(end-5));
                    
                    if (ww~=1) && (interval*(1/Fs)<2)
                        HEAD.SACD.INTER.Time        {kk}{jj}{end+1,:} = NheadTime;
                        HEAD.SACD.INTER.Pos         {kk}{jj}{end+1,:} = NheadPos;
                        HEAD.SACD.INTER.Vel         {kk}{jj}{end+1,:} = headVel;
                        
                        PAT.SACD.INTER.Pos          {kk}{jj}{end+1,:} = NpatPos;
                     	PAT.SACD.INTER.Vel          {kk}{jj}{end+1,:} = patVel;
                        
                        HEAD.SACD.INTER.PosErr      {kk}{jj}{end+1,:} = NerrPos;
                        HEAD.SACD.INTER.VelErr      {kk}{jj}{end+1,:} = NerrVel;
                        
                        HEAD.SACD.INTER.IntPosErr   {kk}{jj}{end+1,:} = NerrPos;
                        HEAD.SACD.INTER.IntVelErr   {kk}{jj}{end+1,:} = NerrVel;
                        
                        % INTER-SACCADE TIME
                        Sacs(mm,21) = 1000*interval*(1/Fs);
                        % ERROR INTER-SACCADE: POSITION,VELOCIY
                        Sacs(mm,17) = abs(mean(NerrPos(end-5)));
                        Sacs(mm,18) = mean(errVel(end-5));
                        % INTEGRATED ERROR INTER-SACCADE: POSITION,VELOCIY
                        Sacs(mm,19) = abs(mean(INTerrPos(end-5)));
                        Sacs(mm,20) = mean(INTerrVel(end-5));
                    else
                        % INTER-SACCADE TIME
                        Sacs(mm,17) = nan;
                        Sacs(mm,18) = nan;
                        Sacs(mm,19) = nan;
                        Sacs(mm,20) = nan;
                        Sacs(mm,21) = nan;
                    end
                    
                    % RATE
                    Sacs(mm,22) = length(pks);
                  	Sacs(mm,23) = VelData(pIdx(ww));  	% peak velocity
                    Sacs(mm,24) = VelData(pIdx(ww))>0; 
                    
                    mm = mm + 1;
                end
            else
                Sacs(1,1:21) = nan;
            end
            
            HEAD.SACD.MAIN = [HEAD.SACD.MAIN ; Sacs];

            if showplot.Sacd
            figure (34) % Position
                subplot(rows,ceil(nTrial/rows),pp) ; hold on ; title(['Fly ' num2str(uFly(kk)) ' Vel ' num2str(uVel(jj))])
                    plot( TimeData , PosData , 'k')  % position data
                    plot(Sacs(:,4),Sacs(:,8), 'g*')  % start
                    plot(Sacs(:,5),Sacs(:,9), 'b*')  % peak
                    plot(Sacs(:,6),Sacs(:,10), 'r*') % end
                    box on
                    xlim([0 10])
                    hold off
                    
            figure (35) % Velocity
                subplot(rows,ceil(nTrial/rows),pp) ; hold on ; title(['Fly ' num2str(uFly(kk)) ' Vel ' num2str(uVel(jj))])
                    plot( TimeData , AbsVelData , 'k') % velocity data
                    plot( TimeData , SacdThresh*ones(length(TimeData),1) ,'c--','LineWidth',3) % threshold line
                    plot(Sacs(:,4),Sacs(:,12), 'g*')  % start
                    plot(Sacs(:,5),Sacs(:,13), 'b*')  % peak
                    plot(Sacs(:,6),Sacs(:,14), 'r*')  % end
                    box on
                    xlim([0 10])
                    hold off 
            end  
            pp = pp + 1; % subplot index
        end
    end
end
disp('DONE')

varNames = {'Fly','Vel','Dir','Trial','StartTime','PeakTime','EndTime','Duration','StartPos','PeakPos','EndPos','Amplitude','StartVel','PeakVel',...
    'EndVel','MaxPreVel','PosErr','VelErr','IntPosErr','IntVelErr','InterDuration','Rate','Peak','SaccadeDir'};


HEAD.SACD.TABLE = splitvars(table(HEAD.SACD.MAIN));
HEAD.SACD.TABLE.Properties.VariableNames = varNames;

%% Saccade Statistics %%
%---------------------------------------------------------------------------------------------------------------------------------
clc
G = [];
pp = [4 5 6];
for jj = 1:length(pp)
    [rr,~] = find( HEAD.SACD.MAIN(:,2) == uVel(pp(jj)) );
    G(rr,1) = jj;
end

SS = [8,12,9,11,14,17,19,18,20,21,22];
YY = {  ['Duration (ms)'],...
        ['Amplitude (' char(176) ')'],...
        ['Trigger Position (' char(176) ')'],...
        ['End Position (' char(176) ')'],...
        ['Peak Velocity (' char(176) '/s)'],...
        ['Error (' char(176) ')'],...
        ['Integrated Error (' char(176) ' \ast s)'],...
      	['Velocity Err (' char(176) '/s)'],...
        ['Integrated Velocity Err (' char(176) ')'],...
        ['Inter-Saccade Duration (ms)'],...
        ['Rate (#/10s)']
                                                        };
CC = {[0 0 0.5],[0 0 0],[0 0.7 0 ],[0.7 0 0],[0.5 0.5 0.5],[0.1 0.5 0.5],[0.5 0.1 0.5],[0.9 0.7 0.7],[0.3 0.3 0.9],...
    [0.6 0.1 0.8],[0 0.8 0.8]};
FF = cell(length(SS),1);
AX = cell(length(SS),1);
for ww = 1:length(SS)
    FF{ww} = figure (400+ww) ; clf ; hold on
    FF{ww}.Position = [-100+(100*ww) , 500-(20*ww) , 700 , 500];
    FF{ww}.Color = 'w';
    AX{ww} = gca;
    axis tight
        bx = boxplot(HEAD.SACD.MAIN(:,SS(ww)),G,'Labels',{num2str(uVel(pp))},'Width',0.5,'Symbol','','Whisker',2);
        ax = gca;
        xlabel(['Stimulus Velocity (' char(176) '/s)'])
        ylabel(YY{ww})
        set(gca,'FontSize',28)
        h = get(bx(5,:),{'XData','YData'});
        for k=1:size(h,1)
           patch(h{k,1},h{k,2},CC{ww});
        end
        set(findobj(gcf,'tag','Median'), 'Color', 'w');
        set(findobj(gcf,'tag','Box'), 'Color', 'k');
        set(findobj(gcf,'tag','Upper Whisker'), 'Color', 'k');
        ax.Children = ax.Children([end 1:end-1]);
end

%% Move Stats to Subplots %%
%---------------------------------------------------------------------------------------------------------------------------------
% Prepare subplots
figure (500) ; clf
% set(gcf,'Position',[400 100 1200 1000])
set(gcf,'Position',[400 100 1200 500])
ss = [6 7];
% ss = [1 2 5 11];
nn = length(ss);
rows = 1;
for ww = 1:nn
    NN(ww) = subplot(rows,ceil(nn/rows),ww);
end
properties = {'XLabel','YLabel','XLim','YLim','XTick','XTickLabel','Box','FontName','FontSize'};
% Paste figures on the subplots
for zz = 1:nn
    copyobj(allchild(get(FF{ss(zz)},'CurrentAxes')),NN(zz));
end
for zz = 1:nn
    for ii = 1:numel(properties)
        prop_str = char(properties(ii));
        set(NN(zz),prop_str,get(AX{ss(zz)},prop_str));
    end
end
%% Saccade Trigger Polar Plot %%
%---------------------------------------------------------------------------------------------------------------------------------
figure (500) ; clf
set(gcf,'Position',[700 300 600 600])
set(gcf,'Color','w')
% polarhistogram(deg2rad(HEAD.SACD.MAIN(:,9)),100,'FaceColor','g','FaceAlpha',.9); hold on
polarhistogram(deg2rad(HEAD.SACD.MAIN(:,11)),100,'FaceColor','r','FaceAlpha',.9); hold on
    legend('Start','End','Location','North')
    set(gca,'FontSize',30)
    grid off
    axis tight
    ax = gca;
    ax.Color = 'k';
    thetalim([-20 20])
    thetaticks([-20 -10 0 10 20])
    thetaticklabels({'20','10','0','-10','-20'})
    d = ax.ThetaDir;
    ax.ThetaZeroLocation = 'top';

%% Inter-Saccade Means %%
%---------------------------------------------------------------------------------------------------------------------------------
for kk = 1:nFly
    for jj = 1:nVel
        if ~isempty(HEAD.SACD.INTER.Time{kk}{jj})
            TimeData = padcat(HEAD.SACD.INTER.Time{kk}{jj}{:});
            PosData = padcat(HEAD.SACD.INTER.Pos{kk}{jj}{:});
            VelData = padcat(HEAD.SACD.INTER.Vel{kk}{jj}{:});
            PatData = padcat(PAT.SACD.INTER.Pos{kk}{jj}{:});
        	PatVelData = padcat(PAT.SACD.INTER.Vel{kk}{jj}{:});
            for ii = 1:size(TimeData,1)
                val = TimeData(ii,:);
                idx = isnan(val);
                percent = mean(idx);
                if (percent>0.1)
                    EI = ii;
                    break
                end
            end
            HEAD.FlyMed.SACD.INTER.Time{kk,1}{jj,1}     = nanmedian(TimeData(1:EI,:),2);
            HEAD.FlyMed.SACD.INTER.Pos{kk,1}{jj,1}      = nanmedian(PosData(1:EI,:),2);
            HEAD.FlyMed.SACD.INTER.Vel{kk,1}{jj,1}   	= nanmedian(VelData(1:EI,:),2);
            PAT.FlyMed.SACD.INTER.Pos{kk,1}{jj,1}       = nanmedian(PatData(1:EI,:),2);
        	PAT.FlyMed.SACD.INTER.Vel{kk,1}{jj,1}       = nanmedian(PatVelData(1:EI,:),2);
        end
    end
end

TT1 = (cat(2,HEAD.FlyMed.SACD.INTER.Time{:}));
TT2 = (cat(2,HEAD.FlyMed.SACD.INTER.Pos{:}));
TT3 = (cat(2,PAT.FlyMed.SACD.INTER.Pos{:}));
TT4 = (cat(2,PAT.FlyMed.SACD.INTER.Vel{:}));
TT5 = (cat(2,HEAD.FlyMed.SACD.INTER.Vel{:}));
for jj = 1:nVel
    TimeData    = padcat(TT1{jj,:});
    PosData     = padcat(TT2{jj,:});
    VelData     = padcat(TT5{jj,:});
    PatData     = padcat(TT3{jj,:});
    PatVelData  = padcat(TT4{jj,:});
    
    for ii = 1:size(TimeData,1)
        val = TimeData(ii,:);
        idx = isnan(val);
        percent = mean(idx);
        if (percent>0.80)
            EI = ii;
            break
        else
            temp = isnan(TimeData);
            [xx,yy] = find(temp==1,1,'first');
            EI = xx;
        end
    end
    
    HEAD.GrandMed.SACD.INTER.Time{jj,1}     = nanmedian(TimeData(1:EI,:),2);
    HEAD.GrandMed.SACD.INTER.Pos{jj,1}      = nanmedian(PosData(1:EI,:),2);
    HEAD.GrandMed.SACD.INTER.Vel{jj,1}  	= nanmedian(VelData(1:EI,:),2);
	PAT.GrandMed.SACD.INTER.Pos{jj,1}       = nanmedian(PatData(1:EI,:),2);
	PAT.GrandMed.SACD.INTER.Vel{jj,1}       = nanmedian(PatVelData(1:EI,:),2);
end

%% Normalized Inter-Saccade Interval Position %
%---------------------------------------------------------------------------------------------------------------------------------
figure (103) ; clf
set(gca,'Color','k')
colorIdx = {'b','g','r','r','g','b'};
for jj = [4,3,5,2,6,1]
    for kk = 1:nFly
        for ii = 1:size(HEAD.SACD.INTER.Pos{kk}{jj},1)
            TimeData = HEAD.SACD.INTER.Time{kk}{jj}{ii}; % get time data
            PosData = HEAD.SACD.INTER.Pos{kk}{jj}{ii};   % get position data
            PatData = PAT.SACD.INTER.Pos{kk}{jj}{ii};    % get pattern data
%             VelData = HEAD.Vel{kk}{jj}(:,ii); % get velocity data
            figure (103) ; hold on ; title('Normalized Inter-Saccade Position') ; box on
                hpos = plot(TimeData , PosData , colorIdx{jj} );
%                 hpat = plot(TimeData , PatData , '--w' , 'LineWidth' , 3 );
%                 uistack(hpat,'top')

                xlabel('Time (s)')
                ylabel('Head Angle (deg)')
                xlim([0 0.6])
                ylim([-22 22])  
        end
    end
end

figure (103) ; hold on
colorIdx = {[0 1 1],[0 0.5 0],[0.7 0 0],[0.7 0 0],[0 0.5 0],[0 1 1]};
for kk = 1:nFly
    for jj = 1:nVel
%         plot(HEAD.FlyMed.SACD.INTER.Time{kk}{jj} , HEAD.FlyMed.SACD.INTER.Pos{kk}{jj} , 'Color', colorIdx{jj} , 'LineWidth',4)
        if kk==nFly
            plot(HEAD.GrandMed.SACD.INTER.Time{jj} , PAT.GrandMed.SACD.INTER.Pos{jj} , '--w' , 'LineWidth',5)
            plot(HEAD.GrandMed.SACD.INTER.Time{jj} , HEAD.GrandMed.SACD.INTER.Pos{jj} , 'Color', 'k' , 'LineWidth',7)
        end
    end
end
% for jj = 1:nVel
%     plot(PAT.Time{1}{jj}(:,1), PAT.Pos{1}{jj}(:,1) , '--w' , 'LineWidth',5)
% end

%% Normalized Inter-Saccade Interval Position Error %
%---------------------------------------------------------------------------------------------------------------------------------
figure (58) ; clf
set(gca,'Color','k')
colorIdx = {'b','g','r','r','g','b'};
for jj = [4,3,5,2,6,1]
    for kk = 1:nFly
        for ii = 1:size(HEAD.SACD.INTER.Pos{kk}{jj},1)
            TimeData = HEAD.SACD.INTER.Time{kk}{jj}{ii};
            PatData  = PAT.SACD.INTER.Pos{kk}{jj}{ii};  
            ErrData  = HEAD.SACD.INTER.PosErr{kk}{jj}{ii}; % get error data
 
            figure (58) ; hold on ; title('Inter-Saccade Interval Position Error') ; box on
                plot(TimeData , ErrData , colorIdx{jj} )
                xlabel('Time (s)')
                ylabel('Head Position Error (deg)')
                xlim([0 0.5])
                ylim([-60 60])  
        end
    end
end

colorIdx = {[0 0.8 1],[0 0.7 0],[0.8 0 0],[0.8 0 0],[0 0.7 0],[0 0.8 1]};
for jj = [4,3,5,2,6,1]
	plot(HEAD.GrandMed.SACD.INTER.Time{jj} , PAT.GrandMed.SACD.INTER.Pos{jj} , 'Color', colorIdx{jj} , 'LineWidth',6 ,...
        'LineStyle','--')
end

%% Inter-Saccade Interval Velocity Step Response %
%---------------------------------------------------------------------------------------------------------------------------------
figure (59) ; clf ; 
% suptitle('Inter-Saccade Interval Velocity')
colorIdx = {'b','g','r','r','g','b'};
pp = 1;
ax = cell(nVel,1);
% subIdx = [4,3,5,2,6,1];
subIdx = [4,5,6];
for jj = subIdx
    for kk = 1:nFly
        for ii = 1:size(HEAD.SACD.INTER.Pos{kk}{jj},1)
            TimeData = HEAD.SACD.INTER.Time{kk}{jj}{ii};
            VelData  = HEAD.SACD.INTER.Vel{kk}{jj}{ii};
            
            figure (59)
            ax{pp} = subplot(1,3,pp) ; hold on ; title(['Stimulus Vel = ' num2str(3.75*uVel(jj)) ' deg/s'])
            box on ; grid on ; set(gca,'Color','k')
            plot(TimeData , VelData , colorIdx{jj} )
            if (pp==nVel || pp==nVel-1) ; xlabel('Time (s)') ; end
%             if mode(pp,2)==1 ; ylabel('Head Velocity (deg/s)') ; end
        end
    end
    pp = pp + 1;
  	plot(HEAD.GrandMed.SACD.INTER.Time{jj} , HEAD.GrandMed.SACD.INTER.Vel{jj} , 'c' , 'LineWidth',5)
	plot(HEAD.GrandMed.SACD.INTER.Time{jj} , PAT.GrandMed.SACD.INTER.Vel{jj} , '-w' , 'LineWidth',5)
    xlabel('Time (s)')
    if subIdx(pp-1)==subIdx(end) ; ylabel('Head (deg/s)') ; end
end
linkaxes([ax{1},ax{2},ax{3},ax{4},ax{5},ax{6}],'xy')
xlim([0 0.24])
ylim([-100 100])
