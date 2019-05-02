classdef IO_Class
    %% IO_Class:
    %   INPUTS:
    %       In      : input  "Fly" class
    %       Out     : output "Fly" class
    
    properties (GetAccess=private) % properties only avaiable to class
        IN                  = []; % input data
        OUT                 = []; % output data
    end
    
    properties (Constant) % constant properties
    end
   
    properties (Dependent) % dependent properties
    end
    
  	properties % all other properties
        BodeFv              = []; % bode plot frequency vector
        BodeGain            = []; % bode plot gain
        BodePhaseDiff       = []; % bode plot phase difference
        
        IOFreq              = []; % frequency vector for specified frequencies

        IOBodeGain          = []; % bode plot gain
        IOBodePhaseDiff     = []; % bode plot phase difference
        
        Coherence           = []; % coherence
     	CoherenceFV       	= []; % coherence frequency vector
        IOCoherence       	= []; % coherence at specified frequencies
        
        CrossCorr           = []; % cross correlation
        TimeLags            = []; % cross correlation time lags
        MaxCC               = []; % max cross correlation
        TimeDiff            = []; % time delay at max cross correlation

        r                   = []; % correlation r-valuets
    	m                   = []; % linear fit slope
        b                   = []; % linear fit y-intercept

    end
    
    methods
        function obj = IO_Class(In,Out)
            if ~isobject(In) || ~isobject(Out)
                error('Input & Output must be "Fly" objects')
            elseif size(In.X,1)~=size(Out.X,1)
                error('1st dimension of objects X-property must match')
            elseif size(In.X,2)~=size(Out.X,2)
                error('2nd dimension of objects X-property must match')
            elseif ~isequal(In.IOFreq,Out.IOFreq)
                error('Input-Output Frequencies much match')
            end
            
            obj.IN  = In;
            obj.OUT = Out;
            
            obj.BodeFv              = In.Fv;
            obj.BodeGain            = Out.Mag ./ In.Mag;
            obj.BodePhaseDiff       = -(In.Phase - Out.Phase);
          	
            obj.IOFreq              = In.IOFreq;
            obj.IOBodeGain          = Out.IOMag ./ In.IOMag;
            obj.IOBodePhaseDiff     = -(In.IOPhase - Out.IOPhase);
            
            for kk = 1:size(In.X,2)
                [obj.Coherence(:,kk),obj.CoherenceFV(:,kk)] = mscohere(In.X(:,kk) , Out.X(:,kk) ,[],[] , In.Fv(:,kk) , In.Fs);
                [obj.IOCoherence(:,kk)] = Get_IO_Cohr(obj.CoherenceFV(:,kk),obj.Coherence(:,kk),In.IOFreq);

                [obj.CrossCorr(:,kk), obj.TimeLags(:,kk) ,obj.MaxCC(:,kk) ,obj.TimeDiff(:,kk)] ...
                    = CrossCorr(In.X(:,kk),Out.X(:,kk),In.Fs);
                
                [obj.r(1,kk),obj.m(1,kk),obj.b(1,kk)] = regression(In.X(:,kk),Out.X(:,kk),'one');
            end
        end
        
        function [] = PlotIOBode(obj,n,lim)
            % PlotIOBode: plot BODE gain & phase-difference for IO frequencies
            %   INPUTS:
            %       obj     : object instance
            %       n       : derivaties of X to plot, default is all
            %       lim     : x-limit
            if nargin==1
                n   = 1:size(obj.BodeGain,2); % all is default
                lim = max(max(obj.BodeFv));
            elseif nargin<=2
                lim = max(max(obj.Fv));
            end
            
            nn  = length(n);

            figure('Name','BODE') ; clf
            pp = 1;
            for kk = n
                subplot(2,nn,pp) ; hold on ; grid on
                title(['X_' num2str(kk)])
%              	plot(obj.BodeFv,obj.BodeGain(:,kk),'k')
                plot(obj.IOFreq,obj.IOBodeGain(:,kk),'r-o')
                xlim([0 lim])
                ylim([0 1.5])
                
                if pp==1
                    ylabel('Gain')
                end
                
                subplot(2,nn,pp+nn) ; hold on ; grid on
                ylabel('Phase')
%              	plot(obj.BodeFv,obj.BodePhaseDiff(:,kk),'k')
                plot(obj.IOFreq,obj.IOBodePhaseDiff(:,kk),'r-o')
             	xlim([0 lim])

                if pp==(n+1)
                    ylabel('Phase Difference')
                end
                
             	xlabel('Frequency')
                
                pp = pp + 1;
            end            
            hold off
        end 
        
        function [] = PlotBode(obj,n,lim)
            % PlotBode: plot BODE gain & phase-difference
            %   INPUTS:
            %       obj     : object instance
            %       n       : derivaties of X to plot, default is all
            %       lim     : x-limit
            if nargin==1
                n   = 1:size(obj.BodeGain,2); % all is default
                lim = max(max(obj.BodeFv));
            elseif nargin<=2
                lim = max(max(obj.BodeFv));
            end
            
            nn  = length(n);

            figure('Name','BODE') ; clf
            pp = 1;
            for kk = n
                subplot(2,nn,pp) ; hold on ; grid on
                title(['X_' num2str(kk)])
             	plot(obj.BodeFv,obj.BodeGain(:,kk),'k')
                xlim([0 lim])
                ylim([0 1.5])
                
                if pp==1
                    ylabel('Gain')
                end
                
                subplot(2,nn,pp+nn) ; hold on ; grid on
                ylabel('Phase')
             	plot(obj.BodeFv,obj.BodePhaseDiff(:,kk),'k')
             	xlim([0 lim])

                if pp==(n+1)
                    ylabel('Phase Difference')
                end
                
             	xlabel('Frequency')
                
                pp = pp + 1;
            end            
            hold off
        end 
        function [] = PlotCohr(obj,n,lim)
            % PlotCohr: plot coherence (and for IO frequencies)
            %   INPUTS:
            %       obj     : object instance
            %       n       : derivaties of X to plot, default is all
            %       lim     : x-limit
            if nargin==1
                n   = 1:size(obj.BodeGain,2); % all is default
                lim = max(max(obj.BodeFv));
            elseif nargin<=2
                lim = max(max(obj.BodeFv));
            end

            nn  = length(n);

            figure('Name','Coherence') ; clf
            pp = 1;
            for kk = n
                subplot(1,nn,pp) ; hold on ; grid on
                title(['X_' num2str(kk)])
                plot(obj.IOFreq,obj.IOCoherence(:,kk),'-or')
                plot(obj.CoherenceFV,obj.Coherence(:,kk),'k')
                xlim([0 lim])
                ylim([0 1])
                if pp==1
                    ylabel('Coherence')
                end
                xlabel('Frequency')

                pp = pp + 1;
            end            
            hold off
        end
    
        function [] = PlotCross(obj,n)
            % PlotCross: plot cross-correlation
            %   INPUTS:
            %       obj     : object instance
            %       n       : derivaties of X to plot, default is all
            if nargin==1
                n   = 1:size(obj.BodeGain,2); % all is default
            end

            nn  = length(n);

            figure('Name','Cross-Correlation') ; clf
            pp = 1;
            for kk = n
                subplot(1,nn,pp) ; hold on ; grid on
                title(['X_' num2str(kk)])
                plot(obj.TimeLags,obj.CrossCorr(:,kk),'-k')
                plot(obj.TimeDiff(:,kk),obj.MaxCC(:,kk),'or')
                tLabel = sprintf('Time Delay = %1.4d',obj.TimeDiff(:,kk));
                text(obj.TimeDiff(:,kk),obj.MaxCC(:,kk),tLabel,'Color','b')
                if pp==1
                    ylabel('Cross-Corelation')
                end
                xlabel('TimeLags')

                pp = pp + 1;
            end            
            hold off
        end
                
        function [] = PlotCorr(obj,n)
            % PlotCorr: plot correlation sactter plot
            %   INPUTS:
            %       obj     : object instance
            %       n       : derivaties of X to plot, default is all
            if nargin==1
                n   = 1:size(obj.BodeGain,2); % all is default
            end

            nn  = length(n);

            figure('Name','Correlation') ; clf
            pp = 1;
            for kk = n
                subplot(1,nn,pp) ; hold on ; grid on
                title(['X_' num2str(kk)])

                lim = max(abs([obj.IN.X(:,kk);obj.OUT.X(:,kk)]));
                xFit = linspace(-lim,lim,10);
                bFit = obj.b(kk) + obj.m(kk)*xFit;

                h = scatplot(obj.IN.X(:,kk),obj.OUT.X(:,kk));
                delete(h.cb)
                plot(xFit,bFit,'r','LineWidth',2)
                text(xFit(7),1.3*bFit(7),['r = ' num2str(obj.r(kk))])
                if pp==1
                    ylabel('OUT')
                end
                xlabel('IN')
                axis(lim*[-1 1 -1 1])

                pp = pp + 1;
            end            
            hold off
        end
    end
end