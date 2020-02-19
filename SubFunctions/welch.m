function [Hs] = welch(x, t, L, overlap, debug)
%% welch: 
%   INPUT:
%       data        :   data to sort
%   OUTPUT:
%       sortData	:   sorted data
%

x = x(:);
t = t(:);

N   = length(x);        % # data points
Ts  = mean(diff(t));    % sampling time
Fs  = round(1 / Ts);    % sampling frequency
Fn = Fs / 2;        	% nyquist frequency

% T   = round(t(end));    % end time

M   = L * Fs;        	% # of points in each segment
S   = M * overlap;   	% # of points to shift between segments
X = buffer(x,M,S);      % partitioned data
T = buffer(t,M,S);      % partitioned time data
K = size(X,2);          % # of segments
W = hann(M);            % Hanning window for each segment

xw = x .* hann(N); % windowed data
Xw = X .* repmat(W,1,K); % windowed segments

% Calculate DFT for windowed and non-windowed data
y = fft(x,[],1) / N;    % normalized fourier transform of raw data
Y = fft(X,[],1) / M;    % normalized fourier transform of partitioned data

yw = fft(xw,[],1) / N;  % normalized fourier transform of windowed raw data
Yw = fft(Xw,[],1) / M;  % normalized fourier transform of windowed partitioned data

Fv  = (linspace(0, 1, fix(N/2)+1)*Fn)';	% frequency vector
Fvp = (linspace(0, 1, fix(M/2)+1)*Fn)';	% frequency vector for partitioned data

Iv  = 1:length(Fv);
Ivp = 1:length(Fvp);

FREQ.y  = y(Iv);     	% complex frequency domain data of raw data
FREQ.Y  = Y(Ivp,:);    	% complex frequency domain data of partitioned data
FREQ.yw = yw(Iv);     	% complex frequency domain data of windowed raw data
FREQ.Yw = Yw(Ivp,:);  	% complex frequency domain data of windowed partitioned data

FREQ.YY = mean(FREQ.Y,2);  	% averaged complex frequency domain data of partitioned data
FREQ.YYw = mean(FREQ.Yw,2);	% averaged complex frequency domain data of windowed partitioned data

if debug
    fig = figure; clf
    set(fig, 'Color', 'w')
    cmap = jet(K);
    plot(t, x, 'Color', 'k', 'LineWidth', 1.5)
    ax = gobjects(K,1);
    for k = 1:K
        ax(k) = subplot(K,1,k) ; hold on
            ylabel(num2str(k))
            plot(t, x, 'Color', 'k', 'LineWidth', 1) 
            plot(T(:,k), X(:,k), 'Color', cmap(k,:), 'LineWidth', 0.75) 
    end
    
    fig = figure; clf ; clear ax
    set(fig, 'Color', 'w')
    ax(1) = subplot(1,2,1) ; hold on
        plot(Fv,  2*abs(FREQ.y), 'Color', 'k', 'LineWidth', 2)
        plot(Fvp, 2*abs(FREQ.Y),'LineWidth', 1)
        plot(Fvp, 2*abs(FREQ.YY), 'r', 'LineWidth', 2)
    ax(2) = subplot(1,2,2) ; hold on
        plot(Fv,  4*abs(FREQ.yw), 'Color', 'k', 'LineWidth', 2)
        plot(Fvp, 4*abs(FREQ.Yw),'LineWidth', 1)
        plot(Fvp, 4*abs(FREQ.YYw), 'r', 'LineWidth', 2)
        
        linkaxes(ax,'xy')
        xlim([0 12])
end


end