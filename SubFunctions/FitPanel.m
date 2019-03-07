function [fitData] = FitPanel(X,t,t_new)
%% FitPanel: fits
%   INPUTS:
%       X       :   discrete arena position data
%       t       :   discrete arena time data
%       t_new  	:   new time vector

%   OUTPUTS:
%       fitData     :   fit data  
%---------------------------------------------------------------------------------------------------------------------------------
% Example Inputs %    
% clear ; clc ; close all    
% data = [];
% t_new = 0:(1/1000):20;
% load('E:\Experiment_HeadExcitation\Chirp\HeadFree\fly_2_trial_8_Amp_15.mat','data','t_p')
% t = t_p;
% X = panel2deg(data(:,2));  % pattern x-pos: subtract mean and convert to deg [deg]
% clear data t_p
%---------------------------------------------------------------------------------------------------------------------------------
    thresh = 3; % velcoity threshold to detect panel transitions 
    n = length(X); % length of signal
    IV = (1:n)'; % index vector
    dX = [abs(diff(X)) ; 0]; % pattern velocity
    mIdx = [IV(1) ; IV(dX>thresh) ; IV(end)]; % detect transition points
    if any(mIdx==n) && any(mIdx==1)
        mIdx = mIdx(2:end-1); % make sure last & 1st point are not detected, we will add it later
    elseif any(mIdx==n)
        mIdx = mIdx(1:end-1); % make sure last point is not detected, we will add it later
	elseif any(mIdx==1)
        mIdx = mIdx(2:end); % make sure first point is not detected, we will add it later 
    end
    mX = X(mIdx); % pattern position corresponding to transtion points
	mt = t(mIdx); % pattern time corresponding to transtion points
    m2Idx = [IV(1) ; mean([mIdx(1:end-1),mIdx(2:end)],2) ; IV(end)]; % peak indicies between transition points
    m2X = X(round(m2Idx)); % peak pattern positions between transition points
  	m2t = t(round(m2Idx)); % peak times between transition points

    [xData, yData] = prepareCurveData( m2t, m2X ); % prepare data for fit
%     ft = 'splineinterp'; % set up fittype and options
    ft = 'linearinterp';

    [fitPattern, ~] = fit( xData, yData, ft, 'Normalize', 'on' ); % fit model to data

    fitData = fitPattern(t_new); % apply fit to new time vector
    
%  	figure (1) ; clf ; hold on
%         plot(t,X,'k')
%     	plot(mt,mX,'g*')
% %         plot(m2t,m2X,'r*')    
%         plot(t_new,fitData,'b')
%         xlim([0 t(end)])
%         xlabel('Time (s)')
%         ylabel('Pattern (deg)')
%         legend('Panel','Peak Marker','Fit') 
end