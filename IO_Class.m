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
        x               = struct('Fv', [], 'Gain', [], 'PhaseDiff', []); % BODE x
        d1x         	= struct('Fv', [], 'Gain', [], 'PhaseDiff', []); % BODE d1x
        d2x         	= struct('Fv', [], 'Gain', [], 'PhaseDiff', []); % BODE d2x
    end
    
    methods
        function obj = IO_Class(IN,OUT)
            
            
        end
        
        function obj = calc_all(obj)
            obj.n               = length(obj.time);                         % # of data points

        end
    end
end