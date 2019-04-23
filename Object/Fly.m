classdef Fly
    %% Fly: computes time & frequency domain quantities
    %   INPUTS:
    %       data    : raw daa
    %       time    : time
    %       Fc      : cutoff frequency
    %       tt      : interpolated time (optional)
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
        Time            = [];       % time
        X               = [];       % filtered raw data
        Fv              = [];       % frequency vector [Hz]
        Mag          	= [];       % magnitude at all frequencies
        Phase        	= [];       % phase at  all frequencies
        IOFreq          = [];       % specified frequencies to find magnitude & gain
        IOMag           = [];       % magnitude at all frequencies
        IOPhase       	= [];       % phase at all frequencies
        Mean            = [];       % mean of data
        STD             = []        % std of data
        AbsMean         = [];       % mean of absolute value of data
        AbsSTD          = [];       % STD of absolute value of data

    end
    
    methods
        function obj = Fly(data,time,Fc,varargin)
            obj.Fc = Fc; % cutoff frequency
            obj.X = data(:); % store data
            obj.Time = time(:); % store time
            
            if nargin==4 % interpolate if new time vector is input
                tt = varargin{1};
                obj.X = interp1(obj.Time, obj.X , tt, 'nearest'); % interpolate to match new time
                obj.Time = tt;
                if min(tt)<min(obj.Time) || max(tt)<max(obj.Time)
                   error('Interpolation time outside of range')
                end
            elseif nargin>4
                error('Too many inputs')
            elseif nargin<2
                error('Not enough many inputs')
            end 
            
            obj = Calc_Main(obj); % run defualt calculations
            
        end
        
        function obj = Calc_Main(obj)
            % Calc_Main: defualt calculations
            obj.n      	= length(obj.Time);                             % # of data points
            obj.Ts     	= mean(diff(obj.Time));                         % sampling Time
            obj.Fs   	= 1/obj.Ts;                                     % sampling frequency
            [b,a]      	= butter(2,obj.Fc/(obj.Fs/2),'low');            % 2nd-order low-pass butterworth filter
            obj.X(:,1) 	= filtfilt(b,a,obj.X);                          % filtered data
            obj.X(:,2)	= filtfilt(b,a,[diff(obj.X(:,1))/obj.Ts ; 0]);	% 1st derivative of data
            obj.X(:,3)	= filtfilt(b,a,[diff(obj.X(:,2))/obj.Ts ; 0]);	% 2nd derivative of data
            
            for kk = 1:size(obj.X,2)
                obj.Mean(1,kk)     	= mean(obj.X(:,kk));             	% mean: data & derivatives
            	obj.AbsMean(1,kk)	= mean(abs(obj.X(:,kk)));          	% mean: absolute value of data & derivatives
                obj.STD(1,kk)     	= std(obj.X(:,kk));                	% std: data & derivatives
                obj.AbsSTD(1,kk)   	= std(obj.X(:,kk));               	% std: absolute value of data & derivatives
                
                [obj.Fv(:,kk),obj.Mag(:,kk),obj.Phase(:,kk)] = FFT(obj.Time,obj.X(:,kk)); % transform data into frequency domain
            end
        end
            
        function obj = IO_Freq(obj,IOFreq)
            % IO_Freq: extract frequency domain data at specified frequencies
            obj.IOFreq = IOFreq(:);
            for kk = 1:size(obj.X,2)
                [obj.IOMag(:,kk),obj.IOPhase(:,kk)] = Get_IO_Freq(obj.Fv(:,kk),obj.Mag(:,kk),obj.Phase(:,kk),obj.IOFreq);
            end
        end
        
        function obj = PlotTime(obj,n)
            % PloTime: plot time domain data
            %   INPUTS:
            %       obj     : object instance
            %       n       : derivaties of X to plot, default is all
            
            if nargin==1
                n = size(obj.X,2);
                n = 1:n;
            end
            
            figure ; clf
            pp = 1;
            for kk = n
                subplot(length(n),1,pp) ; hold on ; grid on
                ylabel(['X_' num2str(kk)])
             	plot(obj.Time,obj.X(:,kk),'k')
                
                if kk==n(end)
                    xlabel('Time')
                end
                
                pp = pp + 1;
            end            
            hold off
        end 
        
	function obj = PlotFreq(obj,n,lim)
            % PlotFreq: plot frequency domain data
            %   INPUTS:
            %       obj     : object instance
            %       n       : derivaties of X to plot, default is all
            
            if nargin==1
                n = size(obj.X,2);
                nn = 1:n;
                lim = max(max(obj.Fv));
            elseif nargin>=2
                nn = size(n);
                lim = max(max(obj.Fv));
            end
            
            figure ; clf
            pp = 1;
            for kk = nn
                subplot(2,length(nn),pp) ; hold on ; grid on
                title(['X_' num2str(kk)])
             	plot(obj.Fv,obj.Mag(:,kk),'k')
                xlim([0 lim])
                
                if pp==1
                    ylabel('Magnitude')
                end
                
                subplot(2,length(nn),pp+n) ; hold on ; grid on
                ylabel(['X_' num2str(kk)])
             	plot(obj.Fv,obj.Phase(:,kk),'k')
             	xlim([0 lim])

                if pp==(n+1)
                    ylabel('Phase')
                end
                
                if kk==nn(end)
                    xlabel('Frequency')
                end
                
                pp = pp + 1;
            end            
            hold off
        end 
        
        
    end
end