
clear ; close all ; clc

% root.Pat = 'Q:\Box Sync\Git\Arena\Patterns';
root.Pat = 'C:\Users\boc5244\Documents\GitHub\Arena\Patterns';
root.Fly = 'H:\EXPERIMENTS\Experiment_Sinusoid\11.25';
root.Func = 'C:\Users\boc5244\Documents\GitHub\Arena\Functions';

% % Pattern
% [file.Pat, path.Pat] = uigetfile({'*.mat', 'MAT-files'}, ...
%     'Select pattern file', root.Pat, 'MultiSelect','off');
% 
% load(fullfile(path.Pat,file.Pat),'pattern');

% % Fly head angles (fly_2_trial_2_freq_2.mat)
% [file.Fly, path.Fly] = uigetfile({'*.mat', 'MAT-files'}, ...
%     'Select fly angles file', fullfile(root.Fly,'Vid','Angles'), 'MultiSelect','off');

load(fullfile(path.Fly,file.Fly),'t_v','hAngles');
load(fullfile(root.Fly,file.Fly),'t_p','data');


%%

% Make eye
delta_phi = 4.6*pi/180; % angle between adjacent ommatidia
n_receptor = 48;        % # of receptors
time_constant = 0.035;  % temproal low-pass filter time constant
Eye = EYE(delta_phi,time_constant,n_receptor); % EYE object

% Pattern
[pattern] = MakePattern_SpatFreq(60);

% Function
[func,~,ftime] = MakePosFunction_Vel(30,2,500);
func = wrap_func(func);

% C = conv2(pattern.Pats(:,:,1),Eye.filt);

[EMD_ALL, Pattern_ALL, Eye_ALL] = EMD(Eye,pattern,func,ftime,false);

[EMD_ALL_filt] = HR_sim(Eye_ALL,ftime);

figure (1) ; clf
imagesc(EMD_ALL_filt)
colormap(jet)

%% EMD Simulation %%
% Make eye
delta_phi = 4.6*pi/180; % angle between adjacent ommatidia
n_receptor = 64;        % # of receptors
time_constant = 0.035;  % temproal low-pass filter time constant
Eye = EYE(delta_phi,time_constant,n_receptor); % EYE object

close all
clear FIG

Vel = flipud(30*[1:6]');
SpatFreq = 3.75*[8]';

maxVel = max(Vel);
T = maxVel/100;
allcomb = combvec(Vel',SpatFreq');
allcomb(3,:) = allcomb(1,:)./allcomb(2,:);
[TempFreq,idx] = sort(allcomb(3,:));
allcomb = round(allcomb(:,idx),3)';
ncomb = size(allcomb,1);
clmn = 3;
hSize = ncomb/clmn;
EMD_Response = nan(ncomb,1);
A = 2;
f = 16;

showplot = true;
if showplot
    FIG(1) = figure (1) ; clf
    FIG(1).Name = 'EMD';
    FIG(1).Position = [100 100 1400 hSize*300];
    movegui(FIG,'center')
    ax = axes;
end

pp = 1;
for jj = 1:ncomb
    vel = allcomb(jj,1);
    spatFreq = allcomb(jj,2);
    tempFreq = allcomb(jj,3);
    
    % Pattern
    [pattern] = MakePattern_SpatFreq(spatFreq);

    % Function
    [func,~,ftime] = MakePosFunction_Vel(vel,T,500);
    func = wrap_func(func);
%     func = A + round(func + A*sin(2*pi*f*ftime));
    
    % Simulate EMD
    [EMD_ALL, Pattern_ALL, Eye_ALL] = EMD(Eye,pattern,func,ftime,false);

    [EMD_ALL_filt] = HR_sim(Eye_ALL,ftime);

    EMD_ALL_filt = EMD_ALL_filt(:,1:end-1); % resize array
    
    EMD_Response(jj) = max(EMD_ALL_filt,[],'all');
    
    if showplot
        % EMD
        ax(pp) = subplot(ceil((ncomb)/clmn),clmn,pp);
        imagesc(EMD_ALL_filt)
        ax(pp).Title.String = ['Temp=' num2str(tempFreq) ' Vel=' num2str(vel) ' Spat=' num2str(spatFreq)];
        ax(pp).YTick = [1 length(ftime)];
        ax(pp).YTickLabels = {'0' num2str(T)};
      	ax(pp).XTick = [1 Eye.n_receptor-1];
        colormap(jet(20))
        ax(pp).XLabel.String = 'Ommatidia';
        rotate3d on
    end

    pp = pp + 1;
end
%%
FIG(2) = figure (2) ; clf ; hold on
FIG(2).Color = 'w';
FIG(2).Position = [100 100 700 500];
ax = gca;
% set(ax,'YColor','w','XColor','w')
ax.XLabel.String = 'Temporal Frequency (cycle/s)';
ax.YLabel.String = 'EMD Response';

plot(TempFreq,EMD_Response,'-or','LineWidth',2)

%% EMD Simulation %%

% Make eye
delta_phi = 4.6*pi/180; % angle between adjacent ommatidia
n_receptor = 16;        % # of receptors
time_constant = 0.035;  % temproal low-pass filter time constant
Eye = EYE(delta_phi,time_constant,n_receptor); % EYE object

close all
clear FIG
Vel = flipud(30*(1:4)');
maxVel = max(Vel);
T = maxVel/100;
SpatFreq = 3.75*[4,8]';
allcomb = combvec(Vel',SpatFreq');
allcomb(3,:) = allcomb(1,:)./allcomb(2,:);
[TempFreq,idx] = sort(allcomb(3,:));
allcomb = round(allcomb(:,idx),2)';
ncomb = size(allcomb,1);
clmn = 6;
hSize = 2*ncomb/clmn;

EMD_Response = nan(ncomb,1);

A = 2;
f = 16;

FIG(1) = figure (1) ; clf
FIG(1).Name = 'EMD';
FIG(1).Position = [100 100 1400 hSize*300];
movegui(FIG,'center')
ax = axes;
showplot = true;
pp = 1;
for jj = 1:ncomb
    vel = allcomb(jj,1);
    spatFreq = allcomb(jj,2);
    tempFreq = allcomb(jj,3);
    
    % Pattern
    [pattern] = MakePattern_SpatFreq(spatFreq);

    % Function
    [func,~,ftime] = MakePosFunction_Vel(vel,T,500);
    func = A + round(func + A*sin(2*pi*f*ftime));
    func = wrap_func(func);

    % Simulate EMD
    [EMD_ALL, Pattern_ALL, Eye_ALL] = EMD(Eye,pattern,func,ftime,false);

    [EMD_ALL_filt] = HR_sim(Eye_ALL);

    EMD_ALL_filt = EMD_ALL_filt(:,1:end-1); % resize array
    
    EMD_Response(jj) = max(EMD_ALL,[],'all');
    
    if showplot
        % Space-time plot
        ax(pp) = subplot(ceil(2*ncomb/clmn),clmn,pp);
        imagesc(Pattern_ALL)
        ax(pp).Title.String = ['Temp=' num2str(tempFreq) ' Vel=' num2str(vel) ' Spat=' num2str(spatFreq)];
        ylabel('time (s)')
        xlabel('Pattern position')
        colormap bone
        freezeColors;

        % EMD
        ax(pp+1) = subplot(ceil((2*ncomb)/clmn),clmn,pp+1);
        imagesc(EMD_ALL_filt);
        colormap(jet);
        xlabel('Ommatidia'); 
        title('EMD')
    end

    pp = pp + 2;
end

%% EMD Simulation %%
close all
clear FIG
vel = flipud(50*(1:5)');
spatFreq = 3.75*[0,2,4,6,8,12,16,24,32,48,96];
for kk = 1:length(vel)
    % Pattern
    [pattern] = MakePattern_SpatFreq(90);
    
    % Function
    [func,~,ftime] = MakePosFunction_Vel(vel(kk),5,100);
    func = wrap_func(func);
    
    % Make eye
    delta_phi = 4.6*pi/180; % angle between adjacent ommatidia
    n_receptor = 20;        % # of receptors
    time_constant = 0.035;  % temproal low-pass filter time constant
    Eye = EYE(delta_phi,time_constant,n_receptor); % EYE object
    
    % Simulate EMD
    [EMD_ALL, Pattern_ALL, Eye_ALL] = EMD(Eye,pattern,func,ftime,false);

    [EMD_ALL_filt] = HR_sim(Eye_ALL);

    EMD_ALL_filt = EMD_ALL_filt(:,1:end-1); % resize array

    FIG(kk) = figure (kk) ; clf
    FIG(kk).Name = ['Vel = ' num2str(vel(kk))];
    if kk>1
        FIG(kk).Position = FIG(kk-1).Position + [50 -50 0 0];
    end
    
    % Space-time plot
    subplot(1,2,1)
    imagesc(Pattern_ALL)
    ylabel('time (s)')
    title('Intensity space-time plot');
    xlabel('Pattern position')
    colormap bone
    freezeColors;

    % EMD
    subplot(1,2,2); 
    imagesc(EMD_ALL_filt);
    colormap(jet);
    xlabel('ommatidia'); 
    title('EMD')
end


%%
close all

tt = TimeVector(t_v,false);
tt = tt(1:end-1);

head.pos = interp1(t_v,hAngles,tt);
head.pos = head.pos - mean(head.pos);
stim.pos = panel2deg(data(:,2));
stim.pos = interp1(t_p,stim.pos,tt);

error.time = tt;
error.pos = stim.pos - head.pos;

% figure (1) ; clf
% hold on
% plot(tt,stim.pos);
% plot(tt,head.pos);
% plot(tt,error.pos);

error.pos   = deg2panel(error.pos) + 20;
stim.pos    = deg2panel(stim.pos) + 20;
