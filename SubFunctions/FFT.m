function [Fv, Mag , Phs , FREQ] = FFT(t,x)
%% FFT: Computes fast-fourier-transform
%   INPUTS:
%       t   :  	time vector [s]
%       x   : 	data
%   OUTPUTS:
%       Fv  :  	frequency vector [Hz]
%       Mag : 	magnitude vector
%       Phs : 	phase vector [rad]
%       FREQ:   complex frequency domain data
%---------------------------------------------------------------------------------------------------------------------------------      
Fs = 1/(mean(diff(t)));                 % sampling frequency [Hz]
L = length(t);                          % length of signal
Fn = Fs/2;                           	% nyquist frequency [Hz]
Fv = (linspace(0, 1, fix(L/2)+1)*Fn)';  % frequency vector [Hz]
fts = fft(x)/L;                        	% normalized fourier transform
Iv = 1:length(Fv);                  	% index vector
Mag = 2*abs(fts(Iv));                   % magnitude
Phs = angle(fts(Iv));                   % phase [rad]
FREQ = fts(Iv);                         % complex frequency domain data
% Phs = rad2deg(Phs);                     % phase [deg]
%---------------------------------------------------------------------------------------------------------------------------------
end