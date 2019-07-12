%% EMD Simulation %%

% Make eye
delta_phi = 4.6*pi/180; % angle between adjacent ommatidia
n_receptor = 16;        % # of receptors
time_constant = 0.035;  % temproal low-pass filter time constant
Eye = EYE(delta_phi,time_constant,n_receptor); % EYE object

close all
clear FIG

% Pick spatial frequencies * vecloities to test
% Vel = flipud(30*[0:0.5:12]');
Vel = [37]';
SpatFreq = 3.75*[8]';

% Make grid
maxVel = max(Vel);
T = maxVel/100;
allcomb = combvec(Vel',SpatFreq');
allcomb(3,:) = allcomb(1,:)./allcomb(2,:);
[TempFreq,idx] = sort(allcomb(3,:));
allcomb = round(allcomb(:,idx),3)';
ncomb = size(allcomb,1);
clmn = 1;
hSize = ncomb/clmn;

EMD_Response = nan(ncomb,1);

showplot = true;
if showplot
    FIG(1) = figure (1) ; clf
    FIG(1).Name = 'EMD';
    FIG(1).Position = [100 100 1400 hSize*250];
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

    [EMD_ALL_filt] = HR(Eye_ALL,ftime);

%     EMD_ALL_filt = EMD_ALL_filt(:,1:end-1); % resize array
    
    EMD_Response(jj) = mean(EMD_ALL_filt,'all');
    
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

FIG(2) = figure (2) ; clf ; hold on
FIG(2).Color = 'w';
FIG(2).Position = [100 100 700 500];
ax = gca;
% set(ax,'YColor','w','XColor','w')
ax.XLabel.String = 'Temporal Frequency (cycle/s)';
ax.YLabel.String = 'EMD Response';

plot(TempFreq,EMD_Response,'-or','LineWidth',2)