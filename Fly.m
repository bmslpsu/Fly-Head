classdef Fly
    %% Fly:
    % 
    properties (GetAccess=private) % properties only avaiable to class
    end
    
    properties (Constant) % constant properties
    end
   
    properties (Dependent) % dependent properties
    end
    
  	properties % all other properties
        n               = [];       % # of data points
        Fc              = [];    	% cutt-off frequency for low pass filter (default=20)
        Ts              = [];       % sampling time
        Fs              = [];       % sampling frequency
        time            = [];       % time
        x               = [];       % filtered raw data
        d1x             = [];       % 1st derivative
        d2x             = [];       % 2nd derivative
        Fv              = [];       % frequency vector [Hz]
        Mag          	= [];       % magnitude at all frequencies
        Phase        	= [];       % phase at  all frequencies
        IOFreq          = [];       % specified frequencies to find magnitude & gain
        IOMag           = [];       % magnitude at all frequencies
        IOPhase       	= [];       % phase at all frequencies
        Mean            = struct('x', [], 'd1x', [] , 'd2x', []); % mean of data
        STD             = struct('x', [], 'd1x', [] , 'd2x', []); % std of data

    end
    
    methods
        function obj = Fly(data,time,Fc,varargin)
            obj.Fc = Fc; % cutoff frequency
            obj.x = data(:); % store time
            obj.time = time(:); % store data
            
            if nargin==4 % interpolate if new time vector is input
                tt = varargin{1};
                obj.x = interp1(obj.time, obj.x , tt, 'nearest'); % interpolate to match new time
                obj.time = tt;
                if min(tt)<min(obj.time) || max(tt)<max(obj.time)
                   error('Interpolation time outside of range')   
                end
            elseif nargin>4
                error('Too many inputs')
            elseif nargin<2
                error('Not enough many inputs')
            end 
            
            obj = calc_all(obj); % run calculations
            
        end
        
        function obj = calc_all(obj)
            obj.n               = length(obj.time);                         % # of data points
            obj.Ts              = mean(diff(obj.time));                     % sampling time
            obj.Fs              = 1/obj.Ts;                                 % sampling frequency
            [b,a]               = butter(2,obj.Fc/(obj.Fs/2),'low');        % 2nd-order low-pass butterworth filter
            obj.x               = filtfilt(b,a,obj.x);                      % filtered data
            obj.d1x             = filtfilt(b,a,[diff(obj.x)/obj.Ts ; 0]); 	% 1st derivative of data
            obj.d2x             = filtfilt(b,a,[diff(obj.d1x)/obj.Ts ; 0]);	% 2nd derivative of data
            obj.Mean.x          = mean(obj.x);                              % mean of data
            obj.Mean.d1x        = mean(abs(obj.d1x));                       % mean of 1st derivative
            obj.Mean.d2x        = mean(abs(obj.d2x));                       % mean of 2nd derivative
            obj.STD.x           = std(obj.x);                               % std of data
            obj.STD.d1x         = std(abs(obj.d1x));                        % std of 1st derivative
            obj.STD.d2x         = std(abs(obj.d2x));                        % std of 2nd derivative

            [obj.Fv,obj.Mag,obj.Phase ] = FFT(obj.time,obj.x);              % transform data into frequency domain
            
            % If specififed, get the frequency domain data at certain frequencies
            if ~isempty(obj.IOFreq)
                obj.IOFreq = obj.IOFreq(:);
                [obj.IOMag,obj.IOPhase] = Get_IO_Freq(obj.Fv,obj.Mag,obj.Phase,obj.IOFreq);
            end
        end
    end
end