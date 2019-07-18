clear;close all;clc

spatPeriod = 3.75*[2,4,6,8,12,16,24,32];
vel = 3.75*[3,5,8,9,11,13,15,20,25,30,40,50];
nVel = length(vel);
nPeriod = length(spatPeriod);
EMD = nan(nPeriod,nVel);
tempFreq = nan(nPeriod,nVel);
% ax = axes;
ncomb = nPeriod*nVel;
clmn = 5;
pp = 1;
for jj = 1:nVel
    for kk = 1:nPeriod
        % Pattern
        [pattern] = MakePattern_SpatFreq(spatPeriod(kk));
        
        tempFreq(kk,jj) = vel(jj)./spatPeriod(kk);

        [EMD_ALL, Pattern_ALL, Eye_ALL] = EMDsim_v5_7_5deg_Grnd_new(pattern, vel(jj), 0, 0);

        [EMD_ALL_filt] = HR_sim(Eye_ALL);

        EMD_ALL_filt = EMD_ALL_filt(:,1:end-1);
        
        EMD(kk,jj) = mean(EMD_ALL_filt,'all');
        
%         % Make space-time plot
%         figure
%         imagesc(Pattern_ALL)
%         ylabel('time (s)');
%         title('Intensity space-time plot');
%         xlabel('Pattern position'); colormap bone; freezeColors;
    
%         % EMD
%         figure
%         imagesc(EMD_ALL_filt);
%         colormap(jet);
%         xlabel('ommatidia')
%         title('EMD')
        
        pp =  pp + 1;
    end
end

%%
figure (2); clf ; hold on
for jj = 1:nVel
    plot([tempFreq(:,jj); 0],[EMD(:,jj) ; 0],'-+','LineWidth',1)
    [mm(jj),idx] = max(EMD(:,jj));
    tmp(jj) = tempFreq(idx,jj);
end
%plot([0,tmp],[0,mm],'o-k','LineWidth',2)
xlim([0 15])


% figure ; clf ; hold on
% for kk = 1:nPeriod
%     plot([0 , tempFreq(kk,:)],[0, EMD(kk,:)],'LineWidth',1)
% end
% xlim([0 15])






