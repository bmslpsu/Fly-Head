%% Select File(s) %%
%---------------------------------------------------------------------------------------------------------------------------------
clear;close all;clc

root = 'Q:\Box Sync\Research\';
[D,I,N,U,T,FILES,PATH] = GetFileData(root);

%% Load video data & run tracking software %%
%---------------------------------------------------------------------------------------------------------------------------------
prev = 0;
WBA = nan(N.file,3);
for kk = 1:N.file
    load([PATH FILES{kk}],'vidData')
    vid = squeeze(vidData);
    dim = size(vid);
    med = median(vid,3);
    
    if (I.fly(kk) - prev)~=0
        disp('Making Mask...')
        [Mask] = MakeWingMask(med);
        close all
    end
    prev = I.fly(kk);

    [Wing] = WingTracker_Edge(med, Mask, 0.35, true);
    
    beep on
    for jj = 1:5
        beep; pause (0.5)
    end
    pause

    WBA(kk,1) = Wing.L;
    WBA(kk,2) = Wing.R;
    
end
close all
WBA(:,3) = WBA(:,1) + WBA(:,2);




