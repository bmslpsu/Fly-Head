function [] = WBA()
%% WBA:
%   INPUTS:
%       -
%
%   OUTPUTS:
%       -
%

% Data location
rootdir = 'H:\EXPERIMENTS\RIGID\Experiment_SOS_v2\tracked_head_wing';

% Select files
[D,I,N,U,T,FILES,PATH,basename] = GetFileData(rootdir,'*.csv',false);

%% Get Data %%
disp('Loading...')
% showplot = false;
figure (1) ; clf
ax(1) = subplot(1,1,1); hold on
set(ax,'LineWidth',1)
swba = nan(N.file,2); ylim([80 180]) ; xlim([0 2100])
[b,a] = butter(2,0.1,'low');
for kk = 1:N.file
    benifly = ImportBenifly(fullfile(PATH, FILES{kk})); % load head angles
    n = size(benifly,1);
    
    lwing = 90 + hampel(1:n,rad2deg(benifly.LWing));
    rwing = 90 + hampel(1:n,rad2deg(benifly.RWing));
    lwing = filtfilt(b,a,lwing);
    rwing = filtfilt(b,a,rwing);
    
	med_lwing = median(lwing);
    med_rwing = median(rwing);
    
    cla
    plot(lwing)
    plot(rwing)
    plot(med_lwing*ones(n,1),'g','LineWidth',1.5)
    plot(med_rwing*ones(n,1),'c','LineWidth',1.5)
    

    swba(kk,1) = med_lwing + med_rwing;
    swba(kk,2) = I.fly(kk);
    %pause(0.1)
    
end
%%
bx = boxplot(swba(:,1), swba(:,2));
test = grpstats(swba(:,1),swba(:,2),'median');
med_all = median(test)
end