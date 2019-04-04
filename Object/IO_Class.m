classdef IO_Class
    %% IO_Class:
    % 
    properties (GetAccess=private) % properties only avaiable to class
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
            
            obj.BodeFv              = In.Fv;
            obj.BodeGain            = Out.Mag ./ In.Mag;
            obj.BodePhaseDiff       = In.Phase - Out.Phase;
          	
            obj.IOFreq              = In.IOFreq;
            obj.IOBodeGain          = Out.IOMag ./ In.IOMag;
            obj.IOBodePhaseDiff     = In.IOPhase - Out.IOPhase;
            

            for kk = 1:size(In.X,2)
                [obj.Coherence(:,kk),obj.CoherenceFV(:,kk)] = mscohere(In.X(:,kk) , Out.X(:,kk) ,[],[] , In.Fv(:,kk) , In.Fs);
                [obj.IOCoherence(:,kk)] = Get_IO_Cohr(obj.CoherenceFV(:,kk),obj.Coherence(:,kk),In.IOFreq);

                [obj.CrossCorr(:,kk), obj.TimeLags(:,kk) ,obj.MaxCC(:,kk) ,obj.TimeDiff(:,kk)] ...
                    = CrossCorr(In.X(:,kk),Out.X(:,kk),In.Fs);
                
                [obj.r(1,kk),obj.m(1,kk),obj.b(1,kk)] = regression(In.X(:,kk),Out.X(:,kk),'one');
            end
            
            
        end
    end
end