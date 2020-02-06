classdef Fly
    %% Fly: computes time domain, frequency domain quantities
    %   INPUTS:
    %       data    : raw data
    %       time    : time
    %       Fc      : cutoff frequency
  	%       IOFreq 	: discrete frequency inputs
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
        Med          	= [];       % median of data
        STD             = []        % std of data
        AbsMean         = [];       % mean of absolute value of data
     	AbsMed        	= [];       % median of absolute value of data
        AbsSTD          = [];       % STD of absolute value of data
        SacdThresh     	= [];       % saccade detetcion threshold
    	SacdCount      	= [];       % saccades in data set
        SacdRate      	= [];       % saccade rate
        SACD            = [];       % saccade data table
        fIdx            = [];       % indicies at IO frequencies
        FREQ            = [];       % complex frequency domain data
        IOFREQ        	= [];       % complex frequency domain data at IO frequencies
        WBF             = [];       % 
      	WBA             = [];       % 
        pxx             = [];       % 
    end
    
    methods
        function obj = Fly(data,time,Fc,IOFreq,tt)
            obj.Fc = Fc; % cutoff frequency
            obj.X = data(:); % store data
            obj.Time = time(:); % store time
            
            if (nargin<5) || isempty(tt)
                tt = obj.Time;
            end
            
            if nargin==3
                IOFreq = 0;
            elseif nargin==5 % interpolate if new time vector is input
                if min(tt)<min(obj.Time) || max(tt)>max(obj.Time)
                   warning('Interpolation time outside of range')
                end
                obj.X = interp1(obj.Time, obj.X , tt, 'linear','extrap'); % interpolate to match new time
                obj.Time = tt;
            elseif nargin>5
                error('Too many inputs')
            elseif nargin<3
                error('Not enough many inputs')
            end 
            
            obj = Calc_Main(obj,IOFreq); % run defualt calculations
            
        end
        
        function obj = Calc_Main(obj,IOFreq)           
            % Calc_Main: default calculations
            obj.n      	= length(obj.Time);                             % # of data points
            obj.Ts     	= mean(diff(obj.Time));                         % sampling Time
            obj.Fs   	= 1/obj.Ts;                                     % sampling frequency
            
            if ~isempty(obj.Fc)
            	[b,a]      	= butter(2,obj.Fc/(obj.Fs/2),'low');     	% 2nd-order low-pass butterworth filter
                obj.X(:,1) 	= filtfilt(b,a,obj.X);                   	% filtered data
            end
            obj.X(:,2)	= [diff(obj.X(:,1))/obj.Ts ; 0];                % 1st derivative of data
            obj.X(:,3)	= [diff(obj.X(:,2))/obj.Ts ; 0];                % 2nd derivative of data    
            
            % Stats
            
            [obj.Fv(:,1),obj.Mag(:,1),obj.Phase(:,1), obj.FREQ(:,1)] = ...
                    FFT(obj.Time,obj.X(:,1));
            
            for kk = 1:size(obj.X,2)
                obj.Mean(1,kk)     	= mean(obj.X(:,kk));             	% mean: data & derivatives
            	obj.AbsMean(1,kk)	= mean(abs(obj.X(:,kk)));          	% mean: absolute value of data & derivatives
                obj.Med(1,kk)     	= median(obj.X(:,kk));             	% median: data & derivatives
            	obj.AbsMed(1,kk)	= median(abs(obj.X(:,kk)));        	% median: absolute value of data & derivatives
                obj.STD(1,kk)     	= std(obj.X(:,kk));                	% std: data & derivatives
                obj.AbsSTD(1,kk)   	= std(obj.X(:,kk));               	% std: absolute value of data & derivatives
                
                [~,obj.Mag(:,kk),obj.Phase(:,kk), obj.FREQ(:,kk)] = ...
                    FFT(obj.Time,obj.X(:,kk)); % transform data into frequency domain
            end
            
            % [obj.pxx,~] = pwelch(obj.X(:,1), round(obj.n/2), [], obj.Fv, 100);
            
            % Input-Output frequency data     
            if ~isempty(IOFreq)
                obj = IO_Freq(obj,IOFreq);
            end
            
            % Saccade detetcion        
            [Sacd,thresh,count,rate] = SacdDetect(obj.X(:,1),obj.Time,2.5,false);
            
            obj.SacdThresh	= thresh;
            obj.SacdCount   = count;
            obj.SacdRate    = rate;
            obj.SACD        = mean(table2array(Sacd),1);
            
        end
        
        function [SACD,THRESH,COUNT,RATE,SaccdRmv_X]  = GetSaccade(obj,SD,debug)
            % GetSaccade: get saccade table
            if nargin<3
                debug = false;
                if nargin==1
                    SD = 2.75;
                end
            end
            
            [SACD,THRESH,COUNT,RATE,SaccdRmv_X] = SacdDetect(obj.X(:,1),obj.Time,SD,debug);
        end
        
        function obj = IO_Freq(obj,IOFreq)
            % IO_Freq: extract frequency domain data at discrete frequencies
            obj.IOFreq = IOFreq(:);
          	[obj.IOMag(:,1),obj.IOPhase(:,1),obj.fIdx] = Get_IO_Freq(obj.Fv(:,1),obj.Mag(:,1),obj.Phase(:,1),obj.IOFreq);
            for kk = 2:size(obj.X,2)
                [obj.IOMag(:,kk),obj.IOPhase(:,kk),~] = Get_IO_Freq(obj.Fv(:,1),obj.Mag(:,kk),obj.Phase(:,kk),obj.IOFreq);
            end        
            
            obj.IOFREQ = obj.FREQ(obj.fIdx,:);
        end
        
        function [] = PlotTime(obj,n)
            % PloTime: plot time domain data
            %   INPUTS:
            %       obj     : object instance
            %       n       : derivaties of X to plot, default is all
            
            if nargin==1
                n = 1:size(obj.X,2);
            end
            
            figure('Name','Time Domain') ; clf
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
        
 	function [] = PlotSacd(obj)
            % PloTime: plot time domain data
            %   INPUTS:
            %       obj     : object instance
            %       n       : derivaties of X to plot, default is all
            
            [Sacd,thresh] = SacdDetect(obj.X(:,1),obj.Time,2.5,false);
                       
            figure('Name','Saccade') ; clf
            
                % Position
                subplot(2,1,1) ; hold on ; grid on
                    ylabel('X_1')
                    plot(obj.Time, zeros(obj.n,1),'--','Color',[0.0 0.5 0.5],'LineWidth',1); % upper detection threshold
                    h.pos   = plot(obj.Time,obj.X(:,1),'k');
                    h.start = plot(Sacd{:,7},Sacd{:,10},'g*'); % start
                    h.peak  = plot(Sacd{:,8},Sacd{:,11},'b*'); % peak
                    h.end   = plot(Sacd{:,9},Sacd{:,12},'r*'); % end
                    legend([h.start h.peak h.end],'Start','Peak','End')
                    xlim([obj.Time(1) obj.Time(end)])
             	% Velocity
                subplot(2,1,2) ; hold on ; grid on
                    ylabel('X_2')
                    plot(obj.Time,obj.X(:,2),'k')
                    plot(Sacd{:,7},Sacd{:,13},'g*') % start
                    plot(Sacd{:,8},Sacd{:,14},'b*') % peak
                    plot(Sacd{:,9},Sacd{:,15},'r*') % end
                    plot(obj.Time, thresh*ones(obj.n,1),'m--','LineWidth',2); % upper detection threshold
                    plot(obj.Time,-thresh*ones(obj.n,1),'m--','LineWidth',2); % lower detection threshold
                	plot(obj.Time, zeros(obj.n,1),'--','Color',[0.0 0.5 0.5],'LineWidth',1); % upper detection threshold
                    xlabel('Time')
                    xlim([obj.Time(1) obj.Time(end)])
            hold off
    end
    
	function [] = PlotFreq(obj,n,lim)
            % PlotFreq: plot frequency domain data
            %   INPUTS:
            %       obj     : object instance
            %       n       : derivaties of X to plot, default is all
            %       lim     : x-limit
            if nargin==1
                n   = 1:size(obj.X,2); % all is default
                lim = max(max(obj.Fv));
            elseif nargin<=2
                lim = max(max(obj.Fv));
            end
            
            nn  = length(n);

          	 figure('Name','Frequency Domain') ; clf
            pp = 1;
            for kk = n
                subplot(2,nn,pp) ; hold on ; grid on
                title(['X_' num2str(kk)])
             	plot(obj.Fv,obj.Mag(:,kk),'k')
                plot(obj.IOFreq,obj.IOMag(:,kk),'r*')
                xlim([0 lim])
                
                if pp==1
                    ylabel('Magnitude')
                end
                
                subplot(2,nn,pp+nn) ; hold on ; grid on
                ylabel('Phase')
             	plot(obj.Fv,obj.Phase(:,kk),'k')
                plot(obj.IOFreq,obj.IOPhase(:,kk),'r*')
             	xlim([0 lim])

                if pp==(n+1)
                    ylabel('Phase')
                end
                
             	xlabel('Frequency')
                
                pp = pp + 1;
            end            
            hold off
        end 
        
        
    end
end