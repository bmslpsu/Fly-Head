classdef Fly
    % Fly:
    % 
    properties (GetAccess=private) % properties only avaiable to class
    end
    
    properties (Constant) % constant properties
    end
   
    properties (Dependent) % dependent properties

    end
    
  	properties % all other properties
        n               = [];       % # of data points
        Fc              = 20;    	% cutt-off frequency for low pass filter [Hz]
        Ts              = [];       % sampling time [s]
        Fs              = [];       % sampling frequency [Hz]
        time            = [];       % time [s]
        x               = [];       % filtered raw data
        d1x             = [];       % 1st derivative
        d2x             = [];       % 2nd derivative
        Fv              = [];       % frequency vector [Hz]
        Magnitude       = [];       % magnitude at frequencies
        Phase        	= [];       % phase at frequencies
        Mean            = struct('x', [], 'd1x', [] , 'd2x', []); % mean of data
        STD             = struct('x', [], 'd1x', [] , 'd2x', []); % std of data

    end
    
    methods
        function obj = Fly(data,tt)
          if nargin==2
             if all(isnumeric([data,tt]))
                obj.x = data(:);
                obj.time = tt(:);
             else
                error('Value must be numeric')
             end  
          end
        end
        
        function obj = calc_all(obj)
            obj.n               = length(obj.time);
            obj.Ts              = mean(diff(obj.time));
            obj.Fs              = 1/obj.Ts;
            [b,a]               = butter(2,obj.Fc/(obj.Fs/2),'low'); % 2nd-order low-pass butterworth filter
            obj.x               = filtfilt(b,a,obj.x);
            obj.d1x             = [diff(obj.x)/obj.Ts ; 0];
            obj.d2x             = [diff(obj.d1x)/obj.Ts ; 0];
            obj.Mean.x          = mean(obj.x);
            obj.Mean.d1x        = mean(abs(obj.d1x));
            obj.Mean.d2x        = mean(abs(obj.d2x));
            obj.STD.x           = std(obj.x);
            obj.STD.d1x         = std(abs(obj.d1x));
            obj.STD.d2x         = std(abs(obj.d2x));

            [obj.Fv,obj.Magnitude,obj.Phase ] = FFT(obj.time,obj.x);
        end
    end
end















