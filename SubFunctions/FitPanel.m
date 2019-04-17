function [fitData] = FitPanel(X,t,t_new,varargin)
%% FitPanel: fits
%   INPUTS:
%       X       :   discrete arena position data
%       t       :   discrete arena time data
%       t_new  	:   new time vector
%       debug   :   if true, show plot
%   OUTPUTS:
%       fitData     :   fit data  
%---------------------------------------------------------------------------------------------------------------------------------
Xw = X;
thresh = 3; % velcoity threshold to detect panel transitions 
n = length(Xw); % length of signal
IV = (1:n)'; % index vector
dX = [abs(diff(Xw)) ; 0]; % pattern velocity
mIdx = [IV(1) ; IV(dX>thresh) ; IV(end)]; % detect transition points
if any(mIdx==n) && any(mIdx==1)
    mIdx = mIdx(2:end-1); % make sure last & 1st point are not detected, we will add it later
elseif any(mIdx==n)
    mIdx = mIdx(1:end-1); % make sure last point is not detected, we will add it later
elseif any(mIdx==1)
    mIdx = mIdx(2:end); % make sure first point is not detected, we will add it later 
end
mX = Xw(mIdx); % pattern position corresponding to transtion points
mt = t(mIdx); % pattern time corresponding to transtion points
m2Idx = [IV(1) ; mean([mIdx(1:end-1),mIdx(2:end)],2) ; IV(end)]; % peak indicies between transition points
m2X = Xw(round(m2Idx)); % peak pattern positions between transition points
m2t = t(round(m2Idx)); % peak times between transition points

% Detect transition points not to fit & fix them
dX = diff(m2X);
dX_mean = mean(dX);
dX_std = std(dX);
dX_err_idx = abs(dX)>(dX_mean+2*dX_std) | circshift(abs(dX)>(dX_mean+1*dX_std),1);
dX_err_time = m2t(dX_err_idx);
dX_err = m2X(dX_err_idx);
err_idx = find(dX_err_idx==1);
for kk = 1:length(err_idx)
    rng = 360 - m2X(err_idx(kk));
    if rng>=360/2
        m2X(err_idx(kk)) = min(m2X);
    else
        m2X(err_idx(kk)) = max(m2X);
    end
end

% Fit points
[xData, yData] = prepareCurveData( m2t, m2X ); % prepare data for fit
%     ft = 'splineinterp'; % set up fittype and options
ft = 'linearinterp';

[fitPattern, ~] = fit( xData, yData, ft, 'Normalize', 'on' ); % fit model to data

fitData = fitPattern(t_new); % apply fit to new time vector

% Detect transition points not from fit & fix them
dX = diff(fitData);
dX_mean = mean(dX);
dX_std = std(dX);
dX_err_idx = abs(dX)>(dX_mean+2*dX_std) | circshift(abs(dX)>(dX_mean+2*dX_std),1);
dX_err_time = t_new(dX_err_idx);
dX_err = fitData(dX_err_idx);
err_idx = find(dX_err_idx==1);
for kk = 1:length(err_idx)
	idx = err_idx(kk);
    for jj = kk:length(fitData)
        if ~ismember(jj,err_idx(kk))
            fitData(err_idx(kk)) = fitData(jj);
        else
            break
        end
    end
end

% fitData = rad2deg(wrapToPi(deg2rad(fitData)));

% debug
if nargin==4
    if varargin{1}==true
        figure (1) ; clf ; hold on
            plot(t,Xw,'k')
%                 plot(mt,mX,'g*')
            plot(m2t,m2X,'r*')
            plot(dX_err_time,dX_err,'og')
            plot(t_new,fitData,'b-')
            xlim([0 t(end)])
            xlabel('Time (s)')
            ylabel('Pattern (deg)')
            legend('Panel','Peak Marker','Fit')
    end
end
end